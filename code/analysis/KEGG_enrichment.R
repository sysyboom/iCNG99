# ======== 🚀 KEGG 富集分析脚本（按 ID 过滤 + 强制斜线 + 固定布局）========

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(clusterProfiler)
  library(KEGGREST)
  library(enrichplot)
  library(ggplot2)
  library(grid)
  library(gtable)
})

# --------- 1. 参数设置 ----------
up_excel    <- "/data/new_data/down_genes1_4.xlsx"  
sheet_up    <- 1
gene_col    <- "gene"
organism    <- "cng"
p_cut       <- 0.05
q_cut       <- 0.2
min_genes   <- 3
max_show    <- 15     
output_pdf  <- "/data/TPM/heat/KEGG_down4.pdf"

# --------- 2. 布局尺寸设置 ----------
pdf_w       <- 10     
pdf_h       <- 8      
panel_w     <- 5      

# --------- 3. 数据处理与富集 ----------
up_genes <- read_excel(up_excel, sheet = sheet_up) |>
  pull(!!gene_col) |> as.character() |> na.omit() |> unique()

kegg_keys <- names(keggList(organism)) |> sub(paste0("^", organism, ":"), "", x = _)
up_kegg <- intersect(up_genes, kegg_keys)

kk_up <- enrichKEGG(
  gene          = up_kegg,
  organism      = organism,
  keyType       = "kegg",
  universe      = kegg_keys,
  pvalueCutoff  = p_cut,
  qvalueCutoff  = q_cut
)

# --------- 4. 过滤逻辑 (排除 00999) ----------
overview_maps <- c(
  "map01100", "map01110", "map01120", "map01200", "map01210", 
  "map01212", "map01220", "map01230", "map01240", "map01250", 
  "map01260", "map01270", "map01280"
)

kk_up@result <- kk_up@result %>%
  # 统一 ID 格式为 mapXXXXX
  mutate(PathwayID = sub("^[a-zA-Z]+", "map", ID)) %>%
  # 排除概览通路
  filter(!PathwayID %in% overview_maps) %>%
  # ✅ 直接通过 ID 过滤：排除 map00999
  filter(PathwayID != "map00999") %>%
  filter(Count >= min_genes)

# --------- 5. 绘图与排序 ----------
if (nrow(kk_up@result) == 0) {
  cat("⚠️ 没有显著通路。\n")
} else {
  # (1) 生成基础图
  p <- dotplot(kk_up, showCategory = max_show)
  
  # (2) 强制排序：锁定因子顺序以形成对角斜线
  p$data$Description <- sub(" - .*", "", p$data$Description)
  
  # 按 GeneRatio 从小到大排序数据
  p$data <- p$data[order(p$data$GeneRatio), ]
  
  # 锁定 Description 因子水平，确保 Y 轴顺序对齐 X 轴数值
  p$data$Description <- factor(p$data$Description, levels = unique(p$data$Description))
  
  # (3) 样式设置
  p <- p +
    ggtitle("KEGG Pathway Enrichment") +
    scale_colour_gradient(low = "#F9844A", high = "#43AA8B") +
    scale_size_continuous(name = "Count", range = c(4, 8)) +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      plot.title  = element_text(hjust = 0.5, face = "bold")
    )
  
  # (4) 锁定绘图框宽度
  g <- ggplotGrob(p)
  panels <- grep("panel", g$layout$name)
  panel_index_w <- unique(g$layout$l[panels])
  g$widths[panel_index_w] <- unit(panel_w, "in")
  
  # (5) 输出 PDF
  pdf(output_pdf, width = pdf_w, height = pdf_h)
  grid.newpage()
  grid.draw(g)
  dev.off()
  
  cat("✅ 任务完成！\n- 已直接过滤 map00999\n- 排序已修正为对角斜线\n- 框宽固定为 5 英寸\n")
}