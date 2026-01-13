## ============================================================
## 🚀 最终版：带星号标注 + 固定框高 + 去网格线 + 原始颜色
## ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

## ========= 1. 路径配置 =========
active_file <- "/data/TPM/heat/active_rxns_union_heat.txt"
rxns_file   <- "/data/new/models/paper/merge_after_YPD_heat.xlsx"
flux_file   <- "/data/TPM/heat/heat4_flux/fluxfold_final.xlsx"

## ========= 2. 读入数据与背景构建 =========
active_ids <- readLines(active_file)
rxns_raw <- read_excel(rxns_file, sheet = "RXNS")

exclude_pattern <- "(?i)exchange|transport|sink|accumulation"
rxns_filtered <- rxns_raw %>%
  select(ID, SUBSYSTEM) %>%
  filter(!is.na(ID), ID %in% active_ids) %>%
  mutate(SUBSYSTEM = if_else(is.na(SUBSYSTEM), "Unannotated", SUBSYSTEM)) %>%
  filter(!grepl(exclude_pattern, SUBSYSTEM))

exclude_names <- c(
  "Metabolic pathways", "Biosynthesis of secondary metabolites",
  "Microbial metabolism in diverse environments", "Carbon metabolism",
  "2-Oxocarboxylic acid metabolism", "Fatty acid metabolism",
  "Biosynthesis of amino acids", "Nucleotide metabolism",
  "Biosynthesis of cofactors", "Nitrogen cycle", "Sulfur cycle"
)
pattern_exclude <- paste(exclude_names, collapse = "|")

rxn_subsystems <- rxns_filtered %>%
  filter(!grepl(pattern_exclude, SUBSYSTEM, ignore.case = TRUE)) %>%
  separate_rows(SUBSYSTEM, sep = ";") %>%
  mutate(SUBSYSTEM = trimws(SUBSYSTEM)) %>%
  filter(SUBSYSTEM != "" & SUBSYSTEM != "Unannotated") %>%
  distinct(ID, SUBSYSTEM)

universe_rxns <- unique(rxn_subsystems$ID)
N <- length(universe_rxns)

## ========= 3. 读入差异集合 =========
up_raw   <- read_excel(flux_file, sheet = "up")
down_raw <- read_excel(flux_file, sheet = "down")

up_ids_filt   <- intersect(unique(up_raw$ID),   universe_rxns)
down_ids_filt <- intersect(unique(down_raw$ID), universe_rxns)
K_up <- length(up_ids_filt); K_down <- length(down_ids_filt)

## ========= 4. 统计计算 (先过滤，后校正) =========
subsystem_stats <- rxn_subsystems %>%
  group_by(SUBSYSTEM) %>%
  summarise(M=n_distinct(ID), x_up=n_distinct(ID[ID %in% up_ids_filt]), x_down=n_distinct(ID[ID %in% down_ids_filt]), .groups="drop") %>%
  filter(M >= 3) 

subsystem_stats <- subsystem_stats %>%
  rowwise() %>%
  mutate(
    p_up   = if(K_up > 0) phyper(x_up - 1, m = M, n = N - M, k = K_up, lower.tail = FALSE) else 1,
    p_down = if(K_down > 0) phyper(x_down - 1, m = M, n = N - M, k = K_down, lower.tail = FALSE) else 1
  ) %>%
  ungroup() %>%
  mutate(p_adj_up = p.adjust(p_up, method="BH"), p_adj_down = p.adjust(p_down, method="BH"))

## ========= 5. 排序与显著性标签生成 =========
subsystem_stats <- subsystem_stats %>%
  mutate(enrich_fold_up = if(K_up > 0) (x_up/K_up)/(M/N) else 0,
         enrich_fold_down = if(K_down > 0) (x_down/K_down)/(M/N) else 0)

sig_up <- subsystem_stats %>%
  filter(p_adj_up < 0.05, x_up > 0) %>%
  mutate(direction = "Upregulated", plot_val = enrich_fold_up, padj = p_adj_up) %>%
  arrange(plot_val) 

sig_down <- subsystem_stats %>%
  filter(p_adj_down < 0.05, x_down > 0) %>%
  filter(!(SUBSYSTEM %in% sig_up$SUBSYSTEM & p_adj_down > p_adj_up)) %>%
  mutate(direction = "Downregulated", plot_val = -enrich_fold_down, padj = p_adj_down) %>%
  arrange(plot_val)

combined_df <- bind_rows(sig_down, sig_up) %>%
  mutate(
    stars = case_when(padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""),
    Pathway = factor(SUBSYSTEM, levels = unique(SUBSYSTEM))
  )

## ========= 6. 绘图 =========
rgb_up_color   <- "#F94144" 
  rgb_down_color <- "#277DA1" 
    x_limit <- ifelse(nrow(combined_df) > 0, max(abs(combined_df$plot_val)) * 1.2, 1) # 留出空间给星号
    
    p <- ggplot(combined_df, aes(x = plot_val, y = Pathway, fill = direction)) +
      geom_col(width = 0.8, color = "white", linewidth = 0.1) +
      # --- 添加星号标注 ---
      geom_text(aes(label = stars), 
                hjust = if_else(combined_df$plot_val > 0, -0.2, 1.2), 
                vjust = 0.8, size = 5) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
      scale_fill_manual(values = c("Upregulated" = rgb_up_color, "Downregulated" = rgb_down_color)) +
      scale_x_continuous(limits = c(-x_limit, x_limit), labels = function(x) abs(x)) +
      labs(x = "Fold Enrichment", y = NULL, fill = NULL) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        aspect.ratio = 1.2, 
        legend.position = "bottom"
      )
    print(p)
    
    ## ========= 7. 保存 =========
    if(nrow(combined_df) > 0) {
      ggsave("/data/TPM/heat/heat4_flux/flux_with_stars.pdf", plot = p, width = 8, height = 7)
      cat("✅ 绘图成功！已标注星号并固定框高。\n")
    }