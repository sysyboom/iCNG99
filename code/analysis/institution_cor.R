library(readxl)
library(dplyr)
library(ggplot2)

# ==================== 1. 数据读取与预处理 ====================
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

# 执行 Log2 转换
merged_log <- merged %>%
  mutate(
    log_x = log2(flux_fold_inst1),
    log_y = log2(flux_fold_inst2)
  ) %>%
  filter(!is.infinite(log_x) & !is.infinite(log_y))

# ==================== 2. 统计计算与 P 值格式化 ====================
pearson_res  <- cor.test(merged_log$log_x, merged_log$log_y, method = "pearson")
spearman_res <- cor.test(merged_log$log_x, merged_log$log_y, method = "spearman")

# 动态 P 值判断函数
format_p_value <- function(p) {
  if (p < 2.22e-16) return("< 2.2e-16")
  else return(sprintf("= %.2e", p))
}

p_val_label <- sprintf("Pearson r = %.3f, p %s\nSpearman rho = %.3f, p %s",
                       pearson_res$estimate, format_p_value(pearson_res$p.value),
                       spearman_res$estimate, format_p_value(spearman_res$p.value))

# ==================== 3. 绘图与输出 PDF (正方形 + 整数坐标) ====================
output_file <- "/data/TPM/vivo/flux_corr.pdf"

# 初始化 PDF 设备
pdf(output_file, width = 6, height = 6)

p <- ggplot(merged_log, aes(x = log_x, y = log_y)) +
  # 象限参考线
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray92", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray92", linewidth = 0.5) +
  # 散点：灰色
  geom_point(alpha = 0.35, color = "#B0B0B0", size = 1.5) +
  # 相关线：统一改为绿色 (森林绿)
  # --- 修改处：开启阴影 (se = TRUE) 并设置透明绿色阴影 ---
  geom_smooth(method = "lm", 
              se = TRUE,           # 开启置信区间阴影
              level = 0.95,        # 阴影范围，默认为 95% 置信区间
              color = "#277DA1",   # 回归线颜色
              fill = "#277DA1",    # 阴影填充颜色，设为与线一致
              alpha = 0.15,        # 阴影透明度，建议在 0.1 到 0.2 之间
              linewidth = 1)  +
  # 强制显示整数刻度
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  theme_classic(base_size = 14) + 
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.title = element_text(size = 12, face = "plain"),
    axis.text = element_text(color = "black"),
    plot.title = element_blank(),
    aspect.ratio = 1 # 强制绘图区为 1:1 正方形
  ) +
  labs(
    x = "Dataset1_flux_log2(Fold Change)",
    y = "Dataset2_flux_log2(Fold Change)"
  ) +
  # 统计标注：左对齐
  annotate(
    "text", x = -Inf, y = Inf,
    label = p_val_label,
    hjust = -0.1, vjust = 1.5, size = 4.5, fontface = "plain"
  )

print(p)

# 关闭 PDF 设备
dev.off()

cat("分析完成！PDF 已保存。坐标轴已整数化，形状为正方形。\n")