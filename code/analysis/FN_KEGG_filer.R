library(readxl)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(KEGGREST)
library(stringr)

file_model   <- "/data/new/models/paper/gene_essentiality_comparison.xlsx"
file_gene2ko <- "/data/new_data/gene2ko_kegg_official.tsv"

model <- read_excel(file_model, col_types = "text")

FN_genes <- model$Gene_ID[
  model$model_predicted == "0" & model$experimental_essential == "1"
]
FN_genes <- unique(str_trim(toupper(FN_genes)))
cat(sprintf("FN genes: %d\n", length(FN_genes)))
background_genes <- unique(toupper(str_trim(model$Gene_ID)))
cat(sprintf("background genes: %d\n", length(background_genes)))

gene2ko <- read.delim(file_gene2ko, stringsAsFactors = FALSE)
colnames(gene2ko) <- c("GID", "KO")
gene2ko$GID <- toupper(str_trim(gene2ko$GID))
gene2ko <- dplyr::distinct(gene2ko)

ko2path <- keggLink("pathway", "ko")
ko2path_df <- data.frame(
  Pathway = sub("path:", "", ko2path),
  KO      = sub("ko:", "", names(ko2path)),
  stringsAsFactors = FALSE
)

ko2path_df <- ko2path_df %>% dplyr::filter(grepl("^map", Pathway))

pathway_names <- keggList("pathway")
pathway2name <- data.frame(
  Pathway     = sub("path:", "", names(pathway_names)),
  Description = pathway_names,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(grepl("^map", Pathway))

overview_maps <- c(
  "map01100","map01110","map01120",
  "map01200","map01210","map01212","map01220","map01230","map01240",
  "map01250","map01260","map01270","map01280"
)

block_patterns <- "Metabolic pathways|Biosynthesis of secondary metabolites|Microbial metabolism|Global and overview maps|Carbohydrate metabolism|Amino acid metabolism"

pathway2name <- pathway2name %>%
  dplyr::filter(!Pathway %in% overview_maps) %>%
  dplyr::filter(!grepl(block_patterns, Description)) %>%
  dplyr::filter(!grepl("Diabetic cardiomyopathy", Description, ignore.case = TRUE))
ko2path_df <- ko2path_df %>%
  dplyr::filter(Pathway %in% pathway2name$Pathway)

pathway2gene <- merge(ko2path_df, gene2ko, by = "KO") %>%
  dplyr::select(term = Pathway, gene = GID) %>%
  dplyr::distinct()

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
  stop("no result")
}
cat(sprintf("%d results \n", nrow(ek_FN@result)))

ek_FN@result <- ek_FN@result %>%
  dplyr::filter(!is.na(Description)) %>%
  dplyr::distinct(Description, .keep_all = TRUE)

p <- dotplot(
  ek_FN,
  showCategory = 20,
  color        = "p.adjust",   
  title        = "KEGG Enrichment (FN genes, filtered)"
) +
  theme_bw(base_size = 14) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  scale_color_gradient(
    low  = "#F9844A",  
    high = "#43AA8B",   
    name = "p.adjust"
  )
print(p)

kegg_sig <- as.data.frame(ek_FN@result) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust)
print(utils::head(
  kegg_sig[, c("ID", "Description", "GeneRatio", "p.adjust", "Count")],
  20
))
