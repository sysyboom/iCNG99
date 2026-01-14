import pandas as pd
import numpy as np
import os

names_path = "/mnt/NFS/fengch/TPM/heat/reference_heat.xlsx"  
var_path = "/mnt/NFS/fengch/TPM/heat/heat4_flux/sol.xlsx"  
output_dir = "/mnt/NFS/fengch/TPM/heat/heat4_flux"

threshold_ratio = 0.8   
tolerance = 1e-6       
use_stable_mode = True  

os.makedirs(output_dir, exist_ok=True)

names_df = pd.read_excel(names_path)
var_df = pd.read_excel(var_path)

names_col = [c for c in names_df.columns if 'name' in c.lower() or 'id' in c.lower()][0]
data_cols = [col for col in var_df.columns if col != 'VAR_NAMES']
result_df = pd.DataFrame(columns=['ID'] + [f'{col}_fold' for col in data_cols])

for name in names_df[names_col]:
    name = str(name).strip()
    nf_name = f'NF_{name}'
    perturb_name = f'PERTURB_NF_{name}'

    nf_row = var_df[var_df['VAR_NAMES'].str.strip().str.upper() == nf_name.upper()]
    perturb_row = var_df[var_df['VAR_NAMES'].str.strip().str.upper() == perturb_name.upper()]

    if not nf_row.empty and not perturb_row.empty:
        row_dict = {'ID': name}
        for col in data_cols:
            nf_val = nf_row[col].values[0]
            pert_val = perturb_row[col].values[0]
            if nf_val == 0 and pert_val == 0:
                row_dict[f'{col}_fold'] = np.nan   
            elif nf_val == 0 and pert_val > 0:
                row_dict[f'{col}_fold'] = np.inf   
            elif nf_val > 0 and pert_val == 0:
                row_dict[f'{col}_fold'] = 0.0      
            else:
                row_dict[f'{col}_fold'] = pert_val / nf_val

        result_df = pd.concat([result_df, pd.DataFrame([row_dict])], ignore_index=True)

fold_file = os.path.join(output_dir, "fluxfold_raw.xlsx")
result_df.to_excel(fold_file, index=False)

cols_wo_id = result_df.columns.difference(['ID'])

def classify_direction(row):
    vals = row[cols_wo_id].replace([np.inf, -np.inf], np.nan)

    if use_stable_mode:
        total = len(vals)
        valid_vals = vals[vals > 0].dropna()
    else:
        valid_vals = vals[vals > 0].dropna()
        total = len(valid_vals)

    if total == 0 or len(valid_vals) == 0:
        return None

    ups = np.sum(valid_vals > 1 + tolerance)
    downs = np.sum(valid_vals < 1 - tolerance)
    equals = np.sum(np.isclose(valid_vals, 1.0, atol=tolerance))

    if ups / total >= threshold_ratio:
        return "up"
    elif downs / total >= threshold_ratio:
        return "down"
    elif equals / total >= threshold_ratio:
        return "no"
    else:
        return None

result_df["Direction"] = result_df.apply(classify_direction, axis=1)

df_up = result_df[result_df["Direction"] == "up"].copy()
df_down = result_df[result_df["Direction"] == "down"].copy()
df_no = result_df[result_df["Direction"] == "no"].copy()

def geometric_mean_consistent(row, direction):
    vals = np.array(row.replace([np.inf, -np.inf], np.nan).dropna())
    vals = vals[vals > 0]
    if len(vals) == 0:
        return np.nan

    if direction == "up":
        consistent_vals = vals[vals > 1]
    elif direction == "down":
        consistent_vals = vals[vals < 1]
    elif direction == "no":
        consistent_vals = vals[np.isclose(vals, 1.0, atol=tolerance)]
    else:
        return np.nan

    if len(consistent_vals) == 0:
        return np.nan
    return np.exp(np.mean(np.log(consistent_vals)))

for df, label in [(df_up, "up"), (df_down, "down"), (df_no, "no")]:
    df["flux_fold"] = df[cols_wo_id].apply(lambda r: geometric_mean_consistent(r, label), axis=1)
    df["Valid_Count"] = df[cols_wo_id].apply(lambda r: np.sum(~np.isnan(r.replace([np.inf, -np.inf], np.nan))), axis=1)

final_excel = os.path.join(output_dir, "fluxfold_final.xlsx" if use_stable_mode else "fluxfold_final_sensitive_mode.xlsx")

with pd.ExcelWriter(final_excel, engine="openpyxl") as writer:
    df_up[["ID", "flux_fold", "Valid_Count"]].to_excel(writer, sheet_name="up", index=False)
    df_down[["ID", "flux_fold", "Valid_Count"]].to_excel(writer, sheet_name="down", index=False)
    df_no[["ID", "flux_fold", "Valid_Count"]].to_excel(writer, sheet_name="no", index=False)

