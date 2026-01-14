import pandas as pd

input_excel = "/mnt/NFS/fengch/TPM/drug/drug_TPM1.xlsx"   
output_excel = "/mnt/NFS/fengch/TPM/drug/Gimme_genes.xlsx"  

df = pd.read_excel(input_excel)

gene_col = df.columns[0]
wt_cols = df.columns[1:4]   
mut_cols = df.columns[4:7]  

df["WT_mean"] = df[wt_cols].mean(axis=1)
df["Mut_mean"] = df[mut_cols].mean(axis=1)

df_out = df[[gene_col, "WT_mean", "Mut_mean"]]
df_out.columns = ["gene", "Gal", "Glu"]

df_out.to_excel(output_excel, index=False)
