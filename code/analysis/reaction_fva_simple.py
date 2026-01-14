import re
import os
import pandas as pd
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis

MODEL_CONTROL = "/mnt/NFS/fengch/new/models/paper/merge_after_YPD.xml"
MODEL_TREATED = "/mnt/NFS/fengch/new/models/paper/merge_after_vivo.xml"

UP_XLSX    = "/mnt/NFS/fengch/TPM/vivo/vivo585_up.xlsx"
DOWN_XLSX  = "/mnt/NFS/fengch/TPM/vivo/vivo585_down.xlsx"

OUT_UP_IDX   = "/mnt/NFS/fengch/TPM/vivo/vivo585_up_index.xlsx"
OUT_UP_VAL   = "/mnt/NFS/fengch/TPM/vivo/vivo585_up_value.xlsx"
OUT_DN_IDX   = "/mnt/NFS/fengch/TPM/vivo/vivo585_down_index.xlsx"
OUT_DN_VAL   = "/mnt/NFS/fengch/TPM/vivo/vivo585_down_value.xlsx"

FVA_FRACTION = 0.99    
TOL_ZERO     = 1e-6    
Q_DMID       = 0.7     
Q_DWIDTH     = 0.7     
OVERLAP_MAX  = 0.6     
TOL_NONZERO  = 1e-6     

FVA_ON_SUBSET = True
FVA_PROCESSES = 8      

EXCLUDE_EX_DM_SK = False
def _normalize_fva_cols(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    for c in df.columns:
        lc = c.lower()
        if lc in ("min", "minimum"):
            rename_map[c] = "minimum"
        elif lc in ("max", "maximum"):
            rename_map[c] = "maximum"
    if rename_map:
        df = df.rename(columns=rename_map)
    return df

def main():
    model_control = read_sbml_model(MODEL_CONTROL)
    model_treated = read_sbml_model(MODEL_TREATED)

    def safe_read_excel(path):
        if path is None or not os.path.exists(path):
            print(f"can't find {path} ")
            return pd.DataFrame(columns=["reaction", "value"])
        df = pd.read_excel(path)
        if df.empty:
            print(f"empty: {path}")
            return pd.DataFrame(columns=["reaction", "value"])
        return df

    up_df = safe_read_excel(UP_XLSX)
    down_df = safe_read_excel(DOWN_XLSX)

    for name, df in {"up_df": up_df, "down_df": down_df}.items():
        if not all(c in df.columns for c in ["reaction", "value"]):
            print(f"{name} lack")
            for c in ["reaction", "value"]:
                if c not in df.columns:
                    df[c] = []

    if EXCLUDE_EX_DM_SK:
        pat = re.compile(r"^(EX_|DM_|SK_)")
        up_df   = up_df  [~up_df["reaction"].astype(str).str.match(pat)].copy()
        down_df = down_df[~down_df["reaction"].astype(str).str.match(pat)].copy()

    print(f"intital: up={len(up_df)}, down={len(down_df)}")

    rxn_list = None
    if FVA_ON_SUBSET:
        model_rxn_ids = {r.id for r in model_control.reactions}
        rxn_list = [r for r in set(up_df["reaction"]).union(set(down_df["reaction"])) if r in model_rxn_ids]

    fva_control = flux_variability_analysis(
        model_control,
        reaction_list=rxn_list,
        fraction_of_optimum=FVA_FRACTION,
        processes=FVA_PROCESSES
    )
    fva_treated = flux_variability_analysis(
        model_treated,
        reaction_list=rxn_list,
        fraction_of_optimum=FVA_FRACTION,
        processes=FVA_PROCESSES
    )

    fva_control = _normalize_fva_cols(fva_control)
    fva_treated = _normalize_fva_cols(fva_treated)

    inactive_control = set(fva_control[(fva_control["minimum"].abs() < TOL_ZERO) &
                                       (fva_control["maximum"].abs() < TOL_ZERO)].index)
    inactive_treated = set(fva_treated[(fva_treated["minimum"].abs() < TOL_ZERO) &
                                       (fva_treated["maximum"].abs() < TOL_ZERO)].index)
    inactive_both = inactive_control & inactive_treated

    filtered_up_df   = up_df  [~up_df["reaction"].isin(inactive_both)].copy()
    filtered_down_df = down_df[~down_df["reaction"].isin(inactive_both)].copy()
    print(f"delete non-active: up={len(filtered_up_df)}, down={len(filtered_down_df)}")

    common = fva_control.index.intersection(fva_treated.index)
    loA, hiA = fva_control.loc[common, "minimum"], fva_control.loc[common, "maximum"]
    loB, hiB = fva_treated.loc[common, "minimum"], fva_treated.loc[common, "maximum"]

    dmid   = ((loA + hiA) / 2.0 - (loB + hiB) / 2.0).abs()
    widthA = (hiA - loA).abs()
    widthB = (hiB - loB).abs()
    dwidth = (widthA - widthB).abs()
    
    cut_dmid   = dmid.quantile(Q_DMID)
    cut_dwidth = dwidth.quantile(Q_DWIDTH)

    overlap = (pd.concat([hiA, hiB], axis=1).min(axis=1) -
               pd.concat([loA, loB], axis=1).max(axis=1)).clip(lower=0)
    union = (pd.concat([hiA, hiB], axis=1).max(axis=1) -
             pd.concat([loA, loB], axis=1).min(axis=1)).clip(lower=1e-12)
    keep_overlap = (overlap / union) < OVERLAP_MAX

    farA = pd.concat([loA.abs(), hiA.abs()], axis=1).max(axis=1)
    farB = pd.concat([loB.abs(), hiB.abs()], axis=1).max(axis=1)
    nonzero_either = (farA >= TOL_NONZERO) | (farB >= TOL_NONZERO)

    primary = (dmid >= cut_dmid) | (dwidth >= cut_dwidth)
    keep = primary & keep_overlap & nonzero_either
    keep_rxns = set(keep[keep].index)

    filtered_up_df   = filtered_up_df  [filtered_up_df["reaction"].isin(keep_rxns)].copy()
    filtered_down_df = filtered_down_df[filtered_down_df["reaction"].isin(keep_rxns)].copy()

    print(f"up={len(filtered_up_df)}, down={len(filtered_down_df)}")

    rxn_to_index = {rxn.id: idx for idx, rxn in enumerate(model_control.reactions)}
    filtered_up_df["index"]   = filtered_up_df["reaction"].map(rxn_to_index)
    filtered_down_df["index"] = filtered_down_df["reaction"].map(rxn_to_index)

    filtered_up_df   = filtered_up_df.dropna(subset=["index"])
    filtered_down_df = filtered_down_df.dropna(subset=["index"])

    filtered_up_df["index"]   = filtered_up_df["index"].astype(int) + 1
    filtered_down_df["index"] = filtered_down_df["index"].astype(int) + 1

    filtered_up_df[["index"]].to_excel(OUT_UP_IDX, index=False, header=False)
    filtered_up_df[["value"]].to_excel(OUT_UP_VAL, index=False, header=False)

    filtered_down_df[["index"]].to_excel(OUT_DN_IDX, index=False, header=False)
    filtered_down_df[["value"]].to_excel(OUT_DN_VAL, index=False, header=False)

if __name__ == "__main__":
    main()
