# ---- edm_code bootstrap ----
source_edm_bootstrap <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    bootstrap <- file.path(current, "R", "bootstrap.R")
    if (file.exists(bootstrap)) {
      source(bootstrap)
      return(source_edm_paths(current))
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find edm_code/R/bootstrap.R. Run this script from inside edm_code.", call. = FALSE)
    }
    current <- parent
  }
}
source_edm_bootstrap()
rm(source_edm_bootstrap)
# ----------------------------

pacman::p_load(
  lubridate,
  tidyverse,
  ISOweek,
  rEDM,
  ggplot2,
  ggforce,
  glue,
  stats,
  dplyr,
  gridExtra,
  cowplot,
  rlang,
  # macam,
  macamts,
  rUIC,
  sinaplot,
  ggExtra,
  ggdensity
)
theme_set(theme_cowplot())

df_all <- read.csv("result/FluSub_JP/smap_coef_res/all.csv")
coef_temp <- gather(df_all, key = "coef", value = "value", A_H3N2_cause_A_H1N1,B_cause_A_H1N1,B_cause_A_H3N2)
ggplot(coef_temp, aes(x = coef, y = value)) +
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24)) +
  geom_boxplot() +
  geom_violin(alpha = 0.5) +
  geom_sina(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  # scale_x_discrete(labels = c("AtoB" = "Type A to B", "AH3N2toAH1N1" = "A/H3N2 to A/H1N1", "BtoAH1N1" = "Type B to A/H1N1", "BtoAH3N2" = "Type B to A/H3N2")) +
  labs(x= "Interactions between Influenza Subtypes",y = "Regularized S-map Coefficient")
# ,title = "Effect of Temperature on Influenza Incidence in Subprefectures"

ggsave("result/FluSub_JP/smap_coef_res/all.tiff",
       units = "in", width = 16, height = 7, dpi = 300, compression = 'lzw')


df <- na.omit(df_all$BtoAH3N2)
median(df)
quantile(df, 0.25)
quantile(df, 0.75)