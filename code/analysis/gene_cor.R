# ==================== 0. 加载必要包 ====================
library(DESeq2)
library(readxl)
library(ggplot2)

# ==================== 1. 读取与预处理数据 ====================
expr1_path <- "/data/PAPER/Q.xlsx"
expr2_path <- "/data/new/transcriptome/paper/vivo_rawreads.xlsx"

expr1 <- read_excel(expr1_path)
expr2 <- read_excel(expr2_path)

expr1 <- as.data.frame(expr1)
rownames(expr1) <- expr1[,1]
expr1 <- expr1[,-1]

expr2 <- as.data.frame(expr2)
rownames(expr2) <- expr2[,1]
expr2 <- expr2[,-1]

coldata1 <- data.frame(
  row.names = colnames(expr1),
  condition = factor(c(rep("in_vitro", 5), rep("in_vivo", 5)))
)

coldata2 <- data.frame(
  row.names = colnames(expr2),
  condition = factor(c(rep("control", 3), rep("lung", 3)))
)

# ==================== 2. DESeq2 差异分析 ====================
dds1 <- DESeqDataSetFromMatrix(countData = round(expr1), colData = coldata1, design = ~ condition)
dds2 <- DESeqDataSetFromMatrix(countData = round(expr2), colData = coldata2, design = ~ condition)

dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)

res1 <- results(dds1)
res2 <- results(dds2)

fc_df <- merge(
  data.frame(gene = rownames(res1), log2FC_1 = res1$log2FoldChange),
  data.frame(gene = rownames(res2), log2FC_2 = res2$log2FoldChange),
  by = "gene"
)

fc_df$log2FC_1 <- as.numeric(as.character(fc_df$log2FC_1))
fc_df$log2FC_2 <- as.numeric(as.character(fc_df$log2FC_2))
fc_df <- na.omit(fc_df)

# ==================== 3. 统计计算与动态 P 值处理 ====================
pearson_test  <- cor.test(fc_df$log2FC_1, fc_df$log2FC_2, method = "pearson")
spearman_test <- cor.test(fc_df$log2FC_1, fc_df$log2FC_2, method = "spearman")

format_p_value <- function(p) {
  if (p < 2.22e-16) return("< 2.2e-16")
  else return(sprintf("= %.2e", p))
}

p_val_label <- sprintf("Pearson r = %.3f, p %s\nSpearman rho = %.3f, p %s",
                       pearson_test$estimate, format_p_value(pearson_test$p.value),
                       spearman_test$estimate, format_p_value(spearman_test$p.value))

# ==================== 4. 绘图并输出 PDF (正方形) ====================
output_pdf <- "/data/TPM/vivo/gene_cor.pdf"

# 设置 PDF 为 6x6 英寸，确保物理形状为正方形
pdf(output_pdf, width = 6, height = 6)

p <- ggplot(fc_df, aes(x = log2FC_1, y = log2FC_2)) +
  # 象限参考虚线
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray92", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray92", linewidth = 0.5) +
  # 散点：灰色
  geom_point(alpha = 0.2, color = "#A0A0A0", size = 0.6) +
  # 趋势线：绿色
  # --- 修改处：开启阴影 (se = TRUE) 并设置透明绿色阴影 ---
  geom_smooth(method = "lm", 
              se = TRUE,           # 开启置信区间阴影
              level = 0.95,        # 阴影范围，默认为 95% 置信区间
              color = "#43AA8B",   # 回归线颜色
              fill = "#43AA8B",    # 阴影填充颜色，设为与线一致
              alpha = 0.15,        # 阴影透明度，建议在 0.1 到 0.2 之间
              linewidth = 1) + 
  # 使用默认比例，不再限制原点居中，让数据填充整个画面
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.title = element_text(size = 12, face = "plain"),
    axis.text = element_text(color = "black"),
    plot.title = element_blank(),
    # 确保绘图区长宽比为 1
    aspect.ratio = 1
  ) +
  labs(
    x = "Dataset1 log2(Fold Change)",
    y = "Dataset2 log2(Fold Change)"
  ) +
  # 统计标注：动态定位到左上角 (Inf 代表画面边缘)
  annotate(
    "text", x = -Inf, y = Inf,
    label = p_val_label,
    hjust = -0.1, vjust = 1.5, size = 4.2, fontface = "plain"
  )

print(p)
dev.off()

cat("正方形 PDF 已生成。坐标轴已恢复为数据自适应缩放。\n")