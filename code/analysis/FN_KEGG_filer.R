#!/usr/bin/env Rscript
# =========================================
# 🚀 KEGG 富集分析（FN 基因）去除概览通路 + 去掉人类疾病通路 + 自定义配色
# =========================================

library(readxl)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(KEGGREST)
library(stringr)

# =========================================
# 🧩 1️⃣ 路径设置
# =========================================
file_model   <- "/data/new/models/paper/gene_essentiality_comparison.xlsx"
file_gene2ko <- "/data/new_data/gene2ko_kegg_official.tsv"

# =========================================
# 🧬 2️⃣ 读取模型数据（FN 基因）
# =========================================
model <- read_excel(file_model, col_types = "text")

# FN：实验必需 (1) & 模型预测非必需 (0)
FN_genes <- model$Gene_ID[
  model$model_predicted == "0" & model$experimental_essential == "1"
]
FN_genes <- unique(str_trim(toupper(FN_genes)))
cat(sprintf("✅ 检测到 FN 基因数: %d\n", length(FN_genes)))

# 背景基因：模型中所有基因
background_genes <- unique(toupper(str_trim(model$Gene_ID)))
cat(sprintf("✅ 背景基因数: %d\n", length(background_genes)))

# =========================================
# 🧬 3️⃣ 读取 KEGG 官方 CNAG→KO 映射
# =========================================
gene2ko <- read.delim(file_gene2ko, stringsAsFactors = FALSE)
colnames(gene2ko) <- c("GID", "KO")
gene2ko$GID <- toupper(str_trim(gene2ko$GID))
gene2ko <- dplyr::distinct(gene2ko)

# =========================================
# 🔗 4️⃣ 获取 KEGG Pathway→KO / Pathway 名称
# =========================================
cat("🔍 从 KEGG 获取 pathway→KO 映射...\n")

ko2path <- keggLink("pathway", "ko")
ko2path_df <- data.frame(
  Pathway = sub("path:", "", ko2path),
  KO      = sub("ko:", "", names(ko2path)),
  stringsAsFactors = FALSE
)

# 只保留 mapXXXXX 通路（去掉 koXXXXX）
ko2path_df <- ko2path_df %>% dplyr::filter(grepl("^map", Pathway))

# 获取通路名称
pathway_names <- keggList("pathway")
pathway2name <- data.frame(
  Pathway     = sub("path:", "", names(pathway_names)),
  Description = pathway_names,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(grepl("^map", Pathway))

# =========================================
# 🚫 5️⃣ 去除 KEGG 概览大通路 + 指定疾病通路
# =========================================
overview_maps <- c(
  "map01100","map01110","map01120",
  "map01200","map01210","map01212","map01220","map01230","map01240",
  "map01250","map01260","map01270","map01280"
)

# 大类关键词（概览、全代谢、大通路）
block_patterns <- "Metabolic pathways|Biosynthesis of secondary metabolites|Microbial metabolism|Global and overview maps|Carbohydrate metabolism|Amino acid metabolism"

pathway2name <- pathway2name %>%
  # 去 overview 的 map011xx / map012xx
  dplyr::filter(!Pathway %in% overview_maps) %>%
  # 去掉 overly general 的大类通路
  dplyr::filter(!grepl(block_patterns, Description)) %>%
  # 🔴 再去掉明显与真菌无关的人病通路
  dplyr::filter(!grepl("Diabetic cardiomyopathy", Description, ignore.case = TRUE))

# 只保留名称表中还存在的通路
ko2path_df <- ko2path_df %>%
  dplyr::filter(Pathway %in% pathway2name$Pathway)

cat("✅ 已去除 KEGG 概览及大类通路，并剔除 Diabetic cardiomyopathy 等无关通路。\n")

# =========================================
# 🧮 6️⃣ 构建 pathway→gene 表（term = pathway, gene = CNAG）
# =========================================
pathway2gene <- merge(ko2path_df, gene2ko, by = "KO") %>%
  dplyr::select(term = Pathway, gene = GID) %>%
  dplyr::distinct()

cat(sprintf("✅ 构建 pathway→gene 映射，共 %d 条记录。\n", nrow(pathway2gene)))

# =========================================
# 🧪 7️⃣ 富集分析（FN 基因）
# =========================================
ek_FN <- enricher(
  gene          = FN_genes,
  TERM2GENE     = pathway2gene,
  TERM2NAME     = pathway2name,
  universe      = background_genes,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

if (is.null(ek_FN) || nrow(ek_FN@result) == 0) {
  stop("❌ 没有显著富集结果，请检查输入。")
}

cat(sprintf("✅ 得到 %d 条 KEGG 富集结果（FN 基因）。\n", nrow(ek_FN@result)))

# =========================================
# 🎨 8️⃣ 绘图（自定义配色的 dotplot）
# =========================================
ek_FN@result <- ek_FN@result %>%
  dplyr::filter(!is.na(Description)) %>%
  dplyr::distinct(Description, .keep_all = TRUE)

# 使用 p.adjust 上色：黄(不显著) → 蓝(更显著)
p <- dotplot(
  ek_FN,
  showCategory = 20,
  color        = "p.adjust",   # 以 p.adjust 作为颜色映射
  title        = "KEGG Enrichment (FN genes, filtered)"
) +
  theme_bw(base_size = 14) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  scale_color_gradient(
    low  = "#F9844A",   # p.adjust 较大 → 亮黄
    high = "#43AA8B",   # p.adjust 较小 → 深蓝
    name = "p.adjust"
  )

print(p)

# =========================================
# 📊 9️⃣ 输出显著通路（前20条）
# =========================================
kegg_sig <- as.data.frame(ek_FN@result) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust)

cat("\n📋 FN 基因显著富集通路（前20条）:\n")
print(utils::head(
  kegg_sig[, c("ID", "Description", "GeneRatio", "p.adjust", "Count")],
  20
))
