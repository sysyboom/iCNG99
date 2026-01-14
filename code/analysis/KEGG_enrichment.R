library(readxl)
library(dplyr)
library(clusterProfiler)
library(KEGGREST)
library(enrichplot)
library(ggplot2)
library(grid)
library(gtable)

up_excel    <- "/data/new_data/down_genes1_4.xlsx"  
sheet_up    <- 1
gene_col    <- "gene"
organism    <- "cng"
p_cut       <- 0.05
q_cut       <- 0.2
min_genes   <- 3
max_show    <- 15     
output_pdf  <- "/data/TPM/heat/KEGG_down4.pdf"

pdf_w       <- 10     
pdf_h       <- 8      
panel_w     <- 5      

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

overview_maps <- c(
  "map01100", "map01110", "map01120", "map01200", "map01210", 
  "map01212", "map01220", "map01230", "map01240", "map01250", 
  "map01260", "map01270", "map01280"
)

kk_up@result <- kk_up@result %>%
  mutate(PathwayID = sub("^[a-zA-Z]+", "map", ID)) %>%
  filter(!PathwayID %in% overview_maps) %>%
  filter(PathwayID != "map00999") %>%
  filter(Count >= min_genes)

if (nrow(kk_up@result) == 0) {
  cat("no notable pathway\n")
} else {
  p <- dotplot(kk_up, showCategory = max_show)
  p$data$Description <- sub(" - .*", "", p$data$Description)
  p$data <- p$data[order(p$data$GeneRatio), ]
  p$data$Description <- factor(p$data$Description, levels = unique(p$data$Description))
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

  g <- ggplotGrob(p)
  panels <- grep("panel", g$layout$name)
  panel_index_w <- unique(g$layout$l[panels])
  g$widths[panel_index_w] <- unit(panel_w, "in")

  pdf(output_pdf, width = pdf_w, height = pdf_h)
  grid.newpage()
  grid.draw(g)
  dev.off()
}
