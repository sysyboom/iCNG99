library(readxl)
library(dplyr)
library(ggplot2)

read_and_merge <- function(file_path) {
  sheets <- c("up", "down", "no")
  df_list <- lapply(sheets, function(s) read_excel(file_path, sheet = s))
  
  bind_rows(df_list) %>%
    select(ID, flux_fold) %>%
    mutate(flux_fold = as.numeric(flux_fold)) %>%
    filter(!is.na(flux_fold)) %>%
    group_by(ID) %>% 
    summarise(flux_fold = mean(flux_fold, na.rm = TRUE))
}

file1 <- "/data/TPM/vivo/vivo_flux/fluxfold_final_Q3.xlsx"
file2 <- "/data/TPM/vivo/vivo_flux/fluxfold_final_vivo1.xlsx"

df1 <- read_and_merge(file1)
df2 <- read_and_merge(file2)
merged <- inner_join(df1, df2, by = "ID", suffix = c("_inst1", "_inst2"))

merged_log <- merged %>%
  mutate(
    log_x = log2(flux_fold_inst1),
    log_y = log2(flux_fold_inst2)
  ) %>%
  filter(!is.infinite(log_x) & !is.infinite(log_y))

pearson_res  <- cor.test(merged_log$log_x, merged_log$log_y, method = "pearson")
spearman_res <- cor.test(merged_log$log_x, merged_log$log_y, method = "spearman")

format_p_value <- function(p) {
  if (p < 2.22e-16) return("< 2.2e-16")
  else return(sprintf("= %.2e", p))
}

p_val_label <- sprintf("Pearson r = %.3f, p %s\nSpearman rho = %.3f, p %s",
                       pearson_res$estimate, format_p_value(pearson_res$p.value),
                       spearman_res$estimate, format_p_value(spearman_res$p.value))

output_file <- "/data/TPM/vivo/flux_corr.pdf"

pdf(output_file, width = 6, height = 6)

p <- ggplot(merged_log, aes(x = log_x, y = log_y)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray92", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray92", linewidth = 0.5) +
  geom_point(alpha = 0.35, color = "#B0B0B0", size = 1.5) +
  geom_smooth(method = "lm", 
              se = TRUE,           
              level = 0.95,        
              color = "#277DA1",   
              fill = "#277DA1",    
              alpha = 0.15,       
              linewidth = 1)  +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
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
    x = "Dataset1_flux_log2(Fold Change)",
    y = "Dataset2_flux_log2(Fold Change)"
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = p_val_label,
    hjust = -0.1, vjust = 1.5, size = 4.5, fontface = "plain"
  )

print(p)
dev.off()
