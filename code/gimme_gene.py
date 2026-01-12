import pandas as pd

# ===== 参数设置 =====
input_excel = "/mnt/NFS/fengch/TPM/drug/drug_TPM1.xlsx"   # 输入文件
output_excel = "/mnt/NFS/fengch/TPM/drug/Gimme_genes.xlsx"  # 输出文件

# 读取数据
df = pd.read_excel(input_excel)

# 假设结构是：Gene, WT_rep1, WT_rep2, WT_rep3, Mut_rep1, Mut_rep2, Mut_rep3
gene_col = df.columns[0]
wt_cols = df.columns[1:4]   # WT 三个重复
mut_cols = df.columns[4:7]  # Mut 三个重复

# 计算平均 TPM
df["WT_mean"] = df[wt_cols].mean(axis=1)
df["Mut_mean"] = df[mut_cols].mean(axis=1)

# 构建输出 DataFrame
df_out = df[[gene_col, "WT_mean", "Mut_mean"]]
df_out.columns = ["gene", "Gal", "Glu"]

# 保存
df_out.to_excel(output_excel, index=False)

print(f"✅ 已生成 GIMME 输入文件: {output_excel}")
print(df_out.head())
