library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

df <- read_excel("growth_test.xlsx")

colnames(df) <- c(
  "source",
  "before_refinement",
  "after_refinement_1",
  "after_refinement_2"
)

df <- df %>%
  mutate(
    row_id = row_number(),
    type = ifelse(
      row_id <= 21,
      "Carbon sources (n = 21)",
      "Nitrogen sources (n = 22)"
    ),
    type = factor(
      type,
      levels = c("Carbon sources (n = 21)", "Nitrogen sources (n = 22)")
    )
  )

plot_df <- df %>%
  pivot_longer(
    cols = c(before_refinement, after_refinement_1, after_refinement_2),
    names_to = "stage",
    values_to = "growth"
  ) %>%
  mutate(
    stage = factor(
      stage,
      levels = c(
        "before_refinement",
        "after_refinement_1",
        "after_refinement_2"
      ),
      labels = c(
        "Draft model",
        "Glucose-based\ngap filling",
        "Multi-substrate\ngap filling"
      )
    ),
    growth = factor(growth, levels = c("No", "Yes"))
  )

plot_df <- plot_df %>%
  group_by(type) %>%
  mutate(
    source = factor(source, levels = rev(unique(source)))
  ) %>%
  ungroup()

col_no <- "#D9D9D9"
  col_yes <- "#43AA8B"

  p <- ggplot(plot_df, aes(x = stage, y = source, fill = growth)) +
    geom_point(
      shape = 21,
      size = 3.0,
      color = "grey60",
      stroke = 0.22
    ) +
    facet_grid(
      rows = vars(type),
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    scale_fill_manual(
      values = c("No" = col_no, "Yes" = col_yes),
      name = NULL,
      labels = c("No biomass", "Biomass produced")
    ) +
    scale_x_discrete(position = "top") +
    labs(
      x = NULL,
      y = NULL,
      title = "Carbon / Nitrogen source growth test"
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10.5,
        margin = margin(b = 4)
      ),
      
      axis.text.x = element_text(
        face = "bold",
        size = 7.5,
        color = "grey25",
        margin = margin(b = 2)
      ),
      axis.text.y = element_text(
        size = 6.8,
        color = "black"
      ),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      
      strip.background = element_rect(
        fill = "grey94",
        color = NA
      ),
      strip.text.y.left = element_text(
        face = "bold",
        size = 8,
        angle = 90,
        color = "black"
      ),
      strip.placement = "outside",
      
      legend.position = "bottom",
      legend.text = element_text(size = 7.2),
      legend.key.size = unit(0.35, "cm"),
      legend.margin = margin(t = -2),
      
      panel.spacing.y = unit(0.35, "lines"),
      panel.grid = element_blank(),
      
      plot.margin = margin(4, 4, 4, 4)
    )
  
  print(p)

  ggsave(
    "growth_test_dotplot_grouped.pdf",
    p,
    width = 4.7,
    height = 7.2,
    units = "in",
    device = cairo_pdf
  )
  
