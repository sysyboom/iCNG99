

library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

map_file     <- "/data/CNAG_GO.tsv"                      
go_desc_file <- "/data/GO_basic_Description.txt"        
status_file  <- "/data/new/models/paper/gene_essentiality_comparison.xlsx"
out_prefix   <- "/data/new/FN_from_CNAG_GO_desc"

cnag_go <- readr::read_delim(map_file, delim = "\t", col_types = cols(.default = col_character())) |>
  dplyr::rename_with(~gsub("\\s+", "_", .x)) |>
  dplyr::rename(
    Gene_ID = dplyr::matches("(?i)^gene(_?id)?$"),
    GO_ID   = dplyr::matches("(?i)^go(_?id)?$")
  ) |>
  dplyr::filter(!is.na(Gene_ID), !is.na(GO_ID), grepl("^GO:\\d+$", GO_ID)) |>
  dplyr::distinct(Gene_ID, GO_ID)

go_desc <- readr::read_delim(go_desc_file, delim = "\t", col_types = cols(.default = col_character())) |>
  dplyr::rename_with(~gsub("\\s+", "_", .x)) |>
  dplyr::rename(
    Class       = dplyr::matches("(?i)^class$"),
    GO          = dplyr::matches("(?i)^go(_?ids)?$"),
    Description = dplyr::matches("(?i)^(description|term|name)$")
  ) |>
  dplyr::mutate(Class = toupper(trimws(Class)), GO = trimws(GO)) |>
  dplyr::filter(Class %in% c("BP","MF","CC"), grepl("^GO:\\d+$", GO)) |>
  dplyr::distinct(GO, Class, Description)

status <- readxl::read_excel(status_file)
stopifnot(all(c("Gene_ID","model_predicted","experimental_essential") %in% names(status)))

bg_genes <- intersect(unique(status$Gene_ID), unique(cnag_go$Gene_ID))

fn_genes <- status |>
  dplyr::filter(model_predicted == 0, experimental_essential == 1) |>
  dplyr::pull(Gene_ID) |>
  unique() |>
  intersect(bg_genes)

cat("background genes:", length(bg_genes), "\n")
cat("FN genes:", length(fn_genes), "\n")

make_t2g <- function(onto = c("BP","MF","CC")) {
  onto <- match.arg(onto)
  t2g <- cnag_go |>
    dplyr::filter(Gene_ID %in% bg_genes) |>
    dplyr::inner_join(go_desc |> dplyr::filter(Class == onto),
                      by = c("GO_ID" = "GO")) |>
    dplyr::transmute(GO = GO_ID, Gene = Gene_ID, Description) |>
    dplyr::distinct()
  list(
    term2gene = t2g |> dplyr::select(GO, Gene) |> dplyr::distinct(),
    term2name = t2g |> dplyr::select(GO, Description) |> dplyr::distinct()
  )
}

run_enrich <- function(t2g, genes, bg, pCut = 0.05, qCut = 0.2) {
  if (is.null(t2g) || nrow(t2g$term2gene) == 0L) return(NULL)
  genes_use <- intersect(genes, bg)
  if (length(genes_use) == 0L) return(NULL)
  clusterProfiler::enricher(
    gene          = genes_use,
    universe      = bg,
    TERM2GENE     = t2g$term2gene,
    TERM2NAME     = t2g$term2name,
    pvalueCutoff  = pCut,
    qvalueCutoff  = qCut
  )
}

t2g_bp <- make_t2g("BP"); ego_bp <- run_enrich(t2g_bp, fn_genes, bg_genes)
t2g_mf <- make_t2g("MF"); ego_mf <- run_enrich(t2g_mf, fn_genes, bg_genes)
t2g_cc <- make_t2g("CC"); ego_cc <- run_enrich(t2g_cc, fn_genes, bg_genes)

as_df_with_class <- function(ego, cls) {
  if (is.null(ego)) return(NULL)
  df <- as.data.frame(ego)
  if (!nrow(df)) return(NULL)
  df$Class <- cls
  df
}

df_all <- list(
  as_df_with_class(ego_bp, "BP"),
  as_df_with_class(ego_mf, "MF"),
  as_df_with_class(ego_cc, "CC")
)
df_all <- Filter(Negate(is.null), df_all)
if (length(df_all)) {
  all_res <- dplyr::bind_rows(df_all)
  write.csv(all_res, paste0(out_prefix, "_GO_merged.csv"), row.names = FALSE)
}

if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp))) write.csv(as.data.frame(ego_bp), paste0(out_prefix, "_GO_BP.csv"), row.names = FALSE)
if (!is.null(ego_mf) && nrow(as.data.frame(ego_mf))) write.csv(as.data.frame(ego_mf), paste0(out_prefix, "_GO_MF.csv"), row.names = FALSE)
if (!is.null(ego_cc) && nrow(as.data.frame(ego_cc))) write.csv(as.data.frame(ego_cc), paste0(out_prefix, "_GO_CC.csv"), row.names = FALSE)

plot_dot <- function(ego, title = NULL, topn = 20) {
  if (is.null(ego) || !nrow(as.data.frame(ego))) {
    message("no notable", if (is.null(title)) "" else title)
    return(invisible(NULL))
  }
  the_title <- if (is.null(title)) "" else title
  
  p <- enrichplot::dotplot(ego, showCategory = topn, color = "p.adjust") +
    ggplot2::ggtitle(the_title) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      plot.title  = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::scale_color_gradient(
      low  = "#F9844A", 
      high = "#43AA8B",  
      name = "p.adjust"
    )
  
  print(p)  
}

plot_dot(ego_bp, "GO Biological Process (FN genes)", topn = 20)
plot_dot(ego_mf, "GO Molecular Function (FN)", topn = 20)
#plot_dot(ego_cc, "GO Cellular Component (FN)", topn = 20)

cover_stat <- function(onto) {
  t2g <- make_t2g(onto)
  if (is.null(t2g) || nrow(t2g$term2gene) == 0L) {
    cat(sprintf("[%s] no useful annotation\n", onto)); return(invisible(NULL))
  }
  bg_cov <- length(intersect(bg_genes, unique(t2g$term2gene$Gene)))
  fn_cov <- length(intersect(fn_genes, unique(t2g$term2gene$Gene)))
  cat(sprintf("[%s] banckground genes overlap=%d, FN genes overlap=%d\n", onto, bg_cov, fn_cov))
}
cover_stat("BP"); cover_stat("MF"); cover_stat("CC")
