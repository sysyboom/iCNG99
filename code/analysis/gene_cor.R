library(DESeq2)
library(readxl)
library(ggplot2)

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

pearson_test  <- cor.test(fc_df$log2FC_1, fc_df$log2FC_2, method = "pearson")
spearman_test <- cor.test(fc_df$log2FC_1, fc_df$log2FC_2, method = "spearman")

format_p_value <- function(p) {
  if (p < 2.22e-16) return("< 2.2e-16")
  else return(sprintf("= %.2e", p))
}

p_val_label <- sprintf("Pearson r = %.3f, p %s\nSpearman rho = %.3f, p %s",
                       pearson_test$estimate, format_p_value(pearson_test$p.value),
                       spearman_test$estimate, format_p_value(spearman_test$p.value))

output_pdf <- "/data/TPM/vivo/gene_cor.pdf"

pdf(output_pdf, width = 6, height = 6)

p <- ggplot(fc_df, aes(x = log2FC_1, y = log2FC_2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray92", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray92", linewidth = 0.5) +
  geom_point(alpha = 0.2, color = "#A0A0A0", size = 0.6) +
  geom_smooth(method = "lm", 
              se = TRUE,          
              level = 0.95,       
              color = "#43AA8B",   
              fill = "#43AA8B",    
              alpha = 0.15,       
              linewidth = 1) + 
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.title = element_text(size = 12, face = "plain"),
    axis.text = element_text(color = "black"),
    plot.title = element_blank(),
    aspect.ratio = 1
  ) +
  labs(
    x = "Dataset1 log2(Fold Change)",
    y = "Dataset2 log2(Fold Change)"
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = p_val_label,
    hjust = -0.1, vjust = 1.5, size = 4.2, fontface = "plain"
  )

print(p)
dev.off()
