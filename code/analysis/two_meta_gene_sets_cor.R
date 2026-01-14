library(readxl)
library(DESeq2)
library(ggplot2)

metabolic_genes <- read_excel("/data/PAPER/meta_gene.xlsx")
gene_list <- metabolic_genes$Gene_ID

res1_meta <- res1[rownames(res1) %in% gene_list, ]
res2_meta <- res2[rownames(res2) %in% gene_list, ]

fc_meta <- merge(
  data.frame(gene = rownames(res1_meta), log2FC_1 = res1_meta$log2FoldChange),
  data.frame(gene = rownames(res2_meta), log2FC_2 = res2_meta$log2FoldChange),
  by = "gene"
)

fc_meta$log2FC_1 <- as.numeric(as.character(fc_meta$log2FC_1))
fc_meta$log2FC_2 <- as.numeric(as.character(fc_meta$log2FC_2))
fc_meta <- na.omit(fc_meta)

spearman_test <- cor.test(fc_meta$log2FC_1, fc_meta$log2FC_2, method = "spearman")

format_p_value <- function(p) {
  if (p < 2.22e-16) return("< 2.2e-16")
  else return(sprintf("= %.2e", p))
}

p_val_label <- sprintf("Spearman rho = %.3f, p %s",
                       spearman_test$estimate, format_p_value(spearman_test$p.value))

output_pdf <- "/data/TPM/vivo/meta_gene_cor.pdf"

get_integer_breaks <- function(x) {
  rng <- range(x, na.rm = TRUE)
  seq(floor(rng[1]), ceiling(rng[2]), by = 2) 
}

pdf(output_pdf, width = 6, height = 6)

p <- ggplot(fc_meta, aes(x = log2FC_1, y = log2FC_2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray92", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray92", linewidth = 0.5) +
  geom_point(alpha = 0.4, color = "#B0B0B0", size = 1.2) +
  geom_smooth(method = "lm", 
              se = TRUE,          
              level = 0.95,        
              color = "#43AA8B",   
              fill = "#43AA8B",    
              alpha = 0.15,        
              linewidth = 1) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 6)) + 
  scale_y_continuous(breaks = scales::breaks_pretty(n = 6)) +
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
    x = "Dataset1_meta_log2(Fold Change)",
    y = "Dataset2_meta_log2(Fold Change)"
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = p_val_label,
    hjust = -0.1, vjust = 1.5, size = 4.5, fontface = "plain"
  )

print(p)
dev.off()
