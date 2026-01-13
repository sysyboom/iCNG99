## =========================
## Excel(raw counts) -> TPM
## =========================

# install.packages(c("readxl", "openxlsx", "dplyr"))
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GenomicFeatures"))

library(readxl)
library(openxlsx)
library(dplyr)
library(GenomicFeatures)

## 1) 参数（按你的实际改）
counts_xlsx <- "/data/new/transcriptome/paper/vivo_rawreads.xlsx"  # 你的raw counts Excel
sheet_name  <- 1                                # sheet名或编号
gene_col    <- "Gene"                           # 基因ID列名（如 CNAG_XXXXX）
gtf_file    <- "/data/new_data/genomic.gtf"
out_xlsx    <- "/data/new_data/vivo_TPM.xlsx"

## 2) 读 Excel counts（要求：一列gene + 多列样本count）
df <- read_excel(counts_xlsx, sheet = sheet_name) |> as.data.frame()

stopifnot(gene_col %in% colnames(df))
genes <- df[[gene_col]] |> as.character()

# 样本列：默认除 gene_col 外全是 count（如果你有其它注释列，请在这里手动指定 sample_cols）
sample_cols <- setdiff(colnames(df), gene_col)

count_mat <- df[, sample_cols, drop = FALSE] |>
  mutate(across(everything(), as.numeric)) |>
  as.matrix()

rownames(count_mat) <- genes

## 3) 从GTF计算 gene exon-reduced length（bp）
library(GenomicRanges)
library(IRanges)

txdb <- makeTxDbFromGFF(gtf_file)
exons_by_gene <- exonsBy(txdb, by = "gene")

# 关键：用 GenomicRanges::reduce（不是 dplyr/purrr 的 reduce）
exons_reduced <- GenomicRanges::reduce(exons_by_gene)

# 关键：用 IRanges::width
gene_len_bp <- sum(IRanges::width(exons_reduced))  # exon并集长度(bp)

gene_len_bp <- as.numeric(gene_len_bp)
names(gene_len_bp) <- names(exons_reduced)

## 4) 对齐 gene_id（只保留两边都有的）
common <- intersect(rownames(count_mat), names(gene_len_bp))
count_mat2 <- count_mat[common, , drop = FALSE]
len_kb <- gene_len_bp[common] / 1000

cat("Genes in Excel:", nrow(count_mat), "\n")
cat("Genes with length from GTF:", length(gene_len_bp), "\n")
cat("Genes used for TPM:", length(common), "\n")

## 5) 计算 TPM
# RPK = counts / length(kb)
rpk <- sweep(count_mat2, 1, len_kb, FUN = "/")

# TPM = RPK / sum(RPK) * 1e6 (按列/样本归一化)
tpm <- sweep(rpk, 2, colSums(rpk, na.rm = TRUE), FUN = "/") * 1e6

## 6) 输出 Excel
tpm_df <- data.frame(gene = rownames(tpm), tpm, check.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "TPM")
writeData(wb, "TPM", tpm_df)

# 可选：输出基因长度与对齐统计
addWorksheet(wb, "gene_length_bp")
writeData(wb, "gene_length_bp",
          data.frame(gene = names(gene_len_bp), length_bp = gene_len_bp))

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat("Saved TPM to:", out_xlsx, "\n")

## 7) sanity check：每个样本TPM总和应接近 1e6
print(round(colSums(tpm, na.rm = TRUE), 2))
