import os
import re
import pandas as pd
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis

MODEL_CONTROL = "/mnt/NFS/fengch/new/models/paper/merge_after_YPD.xml"
MODEL_TREATED = "/mnt/NFS/fengch/new/models/paper/merge_after_vivo.xml"

UP_XLSX    = "/mnt/NFS/fengch/TPM/vivo/vivo1_up.xlsx"
DOWN_XLSX  = "/mnt/NFS/fengch/TPM/vivo/vivo1_down.xlsx"

OUT_UP_IDX   = "/mnt/NFS/fengch/TPM/vivo/test_vivo1_up_index.xlsx"
OUT_UP_VAL   = "/mnt/NFS/fengch/TPM/vivo/test_vivo1_up_value.xlsx"
OUT_DN_IDX   = "/mnt/NFS/fengch/TPM/vivo/test_vivo1_down_index.xlsx"
OUT_DN_VAL   = "/mnt/NFS/fengch/TPM/vivo/test_vivo1_down_value.xlsx"

FVA_FRACTION = 0.8
TOL_ZERO     = 1e-6
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

def safe_read_excel(path):
    if path is None or not os.path.exists(path):
        print(f"can't find")
        return pd.DataFrame(columns=["reaction", "value"])
    df = pd.read_excel(path)
    if df.empty:
        print(f"empty")
        return pd.DataFrame(columns=["reaction", "value"])
    for need in ["reaction","value"]:
        if need not in df.columns:
            df[need] = []
    return df[["reaction","value"]].copy()

def main():
    model_control = read_sbml_model(MODEL_CONTROL)
    model_treated = read_sbml_model(MODEL_TREATED)

    up_df   = safe_read_excel(UP_XLSX)
    down_df = safe_read_excel(DOWN_XLSX)

    if EXCLUDE_EX_DM_SK:
        pat = re.compile(r"^(EX_|DM_|SK_)")
        up_df   = up_df  [~up_df["reaction"].astype(str).str.match(pat)].copy()
        down_df = down_df[~down_df["reaction"].astype(str).str.match(pat)].copy()

    print(f"intital:up={len(up_df)}, down={len(down_df)}")

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

    rxn_to_index = {rxn.id: idx for idx, rxn in enumerate(model_control.reactions)}
    filtered_up_df["index"]   = filtered_up_df["reaction"].map(rxn_to_index)
    filtered_down_df["index"] = filtered_down_df["reaction"].map(rxn_to_index)

    filtered_up_df   = filtered_up_df.dropna(subset=["index"])
    filtered_down_df = filtered_down_df.dropna(subset=["index"])

    filtered_up_df["index"]   = filtered_up_df["index"].astype(int) + 1
    filtered_down_df["index"] = filtered_down_df["index"].astype(int) + 1

    if len(filtered_up_df) > 0:
        filtered_up_df[["index"]].to_excel(OUT_UP_IDX, index=False, header=False)
        filtered_up_df[["value"]].to_excel(OUT_UP_VAL, index=False, header=False)
        print(f"selected up: {len(filtered_up_df)} ")
    else:
        print("no up")

    if len(filtered_down_df) > 0:
        filtered_down_df[["index"]].to_excel(OUT_DN_IDX, index=False, header=False)
        filtered_down_df[["value"]].to_excel(OUT_DN_VAL, index=False, header=False)
        print(f"selected down: {len(filtered_down_df)}")
    else:
        print("no down")

if __name__ == "__main__":
    main()
