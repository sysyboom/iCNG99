import pandas as pd
import numpy as np
import os

# ==================== 参数设置 ====================
names_path = "/mnt/NFS/fengch/TPM/heat/reference_heat.xlsx"   # 含 Names 或 ID
var_path = "/mnt/NFS/fengch/TPM/heat/heat4_flux/sol.xlsx"   # 含 VAR_NAMES + 各条件列
output_dir = "/mnt/NFS/fengch/TPM/heat/heat4_flux"

threshold_ratio = 0.8   # 一致性阈值（例如 0.8 表示 80% 条件方向一致）
tolerance = 1e-6        # 数值容差（判断是否≈1）
use_stable_mode = True  # ✅ True = 包含无流条件（稳定模式），False = 仅有效flux（敏感模式）

os.makedirs(output_dir, exist_ok=True)

# ==================== Step 1: 读取数据 ====================
names_df = pd.read_excel(names_path)
var_df = pd.read_excel(var_path)

# 自动识别 Names 或 ID 列
names_col = [c for c in names_df.columns if 'name' in c.lower() or 'id' in c.lower()][0]
print(f"✅ 检测 Names 列: {names_col}")

# 自动识别所有条件列
data_cols = [col for col in var_df.columns if col != 'VAR_NAMES']
print(f"✅ 检测到 {len(data_cols)} 个条件列: {data_cols[:5]} ...")

# 初始化结果 DataFrame
result_df = pd.DataFrame(columns=['ID'] + [f'{col}_fold' for col in data_cols])

# ==================== Step 2: 计算 Fold Change ====================
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

            # 🚨 计算真实 fold（不设 1.0，保留无效为 NaN）
            if nf_val == 0 and pert_val == 0:
                row_dict[f'{col}_fold'] = np.nan   # 无通量，不计入有效
            elif nf_val == 0 and pert_val > 0:
                row_dict[f'{col}_fold'] = np.inf   # 激活
            elif nf_val > 0 and pert_val == 0:
                row_dict[f'{col}_fold'] = 0.0      # 关闭
            else:
                row_dict[f'{col}_fold'] = pert_val / nf_val

        result_df = pd.concat([result_df, pd.DataFrame([row_dict])], ignore_index=True)

# 保存原始 fold change 结果
fold_file = os.path.join(output_dir, "fluxfold_raw.xlsx")
result_df.to_excel(fold_file, index=False)
print(f"✅ 已保存原始 fold change 结果: {fold_file}")

# ==================== Step 3: 按一致性分类（两种模式） ====================
cols_wo_id = result_df.columns.difference(['ID'])

def classify_direction(row):
    vals = row[cols_wo_id].replace([np.inf, -np.inf], np.nan)

    if use_stable_mode:
        # ✅ 稳定模式：所有条件都算入分母（包括 NaN / 0）
        total = len(vals)
        valid_vals = vals[vals > 0].dropna()
    else:
        # ⚡ 敏感模式：只计算有效 flux（>0 且非 NaN）
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

# 应用分类函数
result_df["Direction"] = result_df.apply(classify_direction, axis=1)

# 按分类拆分
df_up = result_df[result_df["Direction"] == "up"].copy()
df_down = result_df[result_df["Direction"] == "down"].copy()
df_no = result_df[result_df["Direction"] == "no"].copy()

print(f"📊 分类结果: 上调 {len(df_up)} 条，下调 {len(df_down)} 条，不变 {len(df_no)} 条")

# ==================== Step 4: 几何平均（仅一致方向部分） ====================
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

# 应用几何平均与有效计数
for df, label in [(df_up, "up"), (df_down, "down"), (df_no, "no")]:
    df["flux_fold"] = df[cols_wo_id].apply(lambda r: geometric_mean_consistent(r, label), axis=1)
    df["Valid_Count"] = df[cols_wo_id].apply(lambda r: np.sum(~np.isnan(r.replace([np.inf, -np.inf], np.nan))), axis=1)

# ==================== Step 5: 输出结果 ====================
final_excel = os.path.join(output_dir, "fluxfold_final.xlsx" if use_stable_mode else "fluxfold_final_sensitive_mode.xlsx")

with pd.ExcelWriter(final_excel, engine="openpyxl") as writer:
    df_up[["ID", "flux_fold", "Valid_Count"]].to_excel(writer, sheet_name="up", index=False)
    df_down[["ID", "flux_fold", "Valid_Count"]].to_excel(writer, sheet_name="down", index=False)
    df_no[["ID", "flux_fold", "Valid_Count"]].to_excel(writer, sheet_name="no", index=False)

print(f"✅ 已保存最终结果文件: {final_excel}")
print(f"📊 模式: {'稳定 (全条件)' if use_stable_mode else '敏感 (仅有效 flux)'}")
