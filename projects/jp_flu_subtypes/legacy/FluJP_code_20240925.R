rm(list = ls())
library(remotes)
# install.packages("devtools")
# remotes::install_github("ha0ye/rEDM")
# remotes::install_github("yutakaos/rUIC")
# remotes::install_github("ong8181/macam")

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
  stringr,
  sinaplot,
  ggExtra,
  ggdensity,
  paletteer
)
# Load library
# packageVersion("macam") # v 0.1.4
packageVersion("rEDM") # v 0.7.5
packageVersion("macamts") # v 0.1.4/ v 0.2.0 2025/08/14
packageVersion("rUIC") # v 0.9.12/ v 0.9.15 2025/08/14
# library(rEDM)
# vignette("rEDM-tutorial")
theme_set(theme_cowplot())


## Set path
setwd("/Users/ayo/Desktop/_GSAIS_/mzmtlab/microbiome dynamics/")
source("code/JP_fluSub/functions_FluSubJP.R")
## result save path
res_save_path <- sprintf("result/FluSub_JP/%s/", Sys.Date())
dir.create(file.path(res_save_path))
print(res_save_path)

## plot save path
fig_date_path <- paste0(res_save_path, "figure/")
dir.create(file.path(fig_date_path))
print(fig_date_path)
dir.create(file.path(paste0(fig_date_path, "smapFig/")))
# Raw time series
raw_ts_path <- paste0(fig_date_path, "raw_ts/")
dir.create(file.path(raw_ts_path))
# Standalized time series
standalized_ts_path <- paste0(fig_date_path, "standarlized_incidence_H1H3B.tiff")
# UIC figure
uic_plot_path <- paste0(res_save_path, "figure/uicFig/")
dir.create(file.path(uic_plot_path))
xmap_plot_path <- paste0(res_save_path, "figure/crossMapping/")
dir.create(file.path(xmap_plot_path))
dir.create(file.path(paste0(uic_plot_path, "uic/")))


# tables
dir.create(file.path(paste0(res_save_path, "table/")))
dir.create(file.path(paste0(res_save_path, "table/smap")))
dir.create(file.path(paste0(res_save_path, "table/smap/block")))
dir.create(file.path(paste0(res_save_path, "table/smap/coef")))
dir.create(file.path(paste0(res_save_path, "table/smap/pred_res")))

uic_surr_tbl_path <- paste0(res_save_path, "table/uic/")
print(uic_surr_tbl_path)
dir.create(file.path(uic_surr_tbl_path))
dir.create(file.path(paste0(uic_surr_tbl_path, "result/")))
dir.create(file.path(paste0(uic_surr_tbl_path, "surrogate_dt/")))
dir.create(file.path(paste0(uic_plot_path, "uic_surr/")))


## --------------------- Preparation --------------------- ##
## import data
# dt0 <- read.csv("data/FluSub_jp/FluSub_10to19_jp.csv")
# dt0 <- dt0[c(70:537),]

dt0 <- read.csv("data/FluSub_jp/FluSub_11to19_jp_per_20240925.csv")

# dt <- subset(dt, select = -c(X))
dt0 <- dt0 %>%
  mutate(date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1")))
dt0$date <- as.Date(dt0$date)

# log(): base "e"
df <- data.frame(Date = dt0$date,
                 A = log(dt0$A + 1),
                 B_Victoria = log(dt0$B_Victoria + 1),
                 B_Yamagata = log(dt0$B_Yamagata + 1),

                 A_H1N1 = log(dt0$A_H1N1 + 1),
                 A_H3N2 = log(dt0$A_H3N2 + 1),
                 B = log(dt0$B + 1)
)

# plot with the log ts
# ggplot(df, aes(Date)) +
ggplot(dt0, aes(date)) +
  geom_point(aes(y = A_H1N1, colour = "A/H1N1")) +
  geom_point(aes(y = A_H3N2, colour = "A/H3N2")) +
  geom_point(aes(y = B, colour = "B")) +

  geom_line(aes(y = A_H1N1, colour = "A/H1N1")) +
  geom_line(aes(y = A_H3N2, colour = "A/H3N2")) +
  geom_line(aes(y = B, colour = "B")) +

  # labs(colour = "Influenza Subtypes",) +
  # scale_colour_grey(start = 0, end = .6) +

  # theme_bw() +
  # ylim(c(0, 10)) +
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24)) +
  xlab("Date") +
  ylab("Incidence")

ggplot(dt0, aes(date)) +
  geom_point(aes(y = A_H1N1, colour = "A/H1N1")) +
  geom_point(aes(y = A_H3N2, colour = "A/H3N2")) +
  geom_point(aes(y = B, colour = "B")) +

  geom_line(aes(y = A_H1N1, colour = "A/H1N1")) +
  geom_line(aes(y = A_H3N2, colour = "A/H3N2")) +
  geom_line(aes(y = B, colour = "B")) +

  labs(colour = "Influenza Subtypes") +
  scale_x_date(date_breaks = "6 month", date_labels = "%Y-%m") + # 设置日期间隔和格式
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1), # 旋转 x 轴标签
        axis.text.y = element_text(size = 24),
        legend.position = "none") +
  scale_color_paletteer_d("fishualize::Acanthurus_sohal") +
  xlab("Date") +
  ylab("Incidence")
ggsave("/Users/ayo/Desktop/_GSAIS_/mzmtlab/microbiome dynamics/result/FluSub_JP/raw_ts.tiff",
       units = "in", width = 38, height = 7, dpi = 300, compression = 'lzw')

# ggsave(paste0("/Users/yining/Desktop/_GSAIS_/2024/学会/グローバルヘルス合同大会/raw_timeseries.tiff"),
#        units = "in", width = 18, height = 8, dpi = 300, compression = 'lzw')
#
# print("---")

## -------------------------------------------------------------------------- ##
## -------------------------- Parameters: E, theta -------------------------- ##
## -------------------------------------------------------------------------- ##
### Simplex Projection & Univariate S-map
### Only output the optimal E&theta result, if want to output the plot,
# see "Flu_subJP_oldCode_paper.R - Simplex Projection, S-map"
# E_range <- 52
# theta_range <- c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
# E_results_df <- data.frame(Variable = character(),
#                            E_value = numeric(),
#                            rho = numeric(),
#                            rmse = numeric(),
#                            stringsAsFactors = FALSE)
# Smap_results_df <- data.frame(Variable = character(),
#                               theta_value = numeric(),
#                               rho = numeric(),
#                               rmse = numeric(),
#                               stringsAsFactors = FALSE)
# for (i in 2:length(df)) {
#   simplex_projection <- rEDM::simplex(df[[colnames(df)[i]]],
#                                       E = 1:E_range,
#                                       silent = T,
#                                       stats_only = FALSE)
#   E_value <- simplex_projection[which.min(simplex_projection$rmse), "E"]
#   simplex_res <- rEDM::simplex(df[[colnames(df)[i]]],
#                                E = E_value,
#                                stats_only = FALSE)
#   E_results_df <- rbind(E_results_df,
#                         data.frame(Variable = colnames(df)[i],
#                                    E_value = simplex_res$E,
#                                    rho = simplex_res$rho,
#                                    rmse = simplex_res$rmse))
#
#   ## S-map
#   uni_Smap <- rEDM::s_map(df[[colnames(df)[i]]],
#                           E = E_value, theta = theta_range,
#                           silent = T,
#                           stats_only = FALSE)
#   theta_value <- uni_Smap[which.min(uni_Smap$rmse), "theta"]
#   smap_res <- rEDM::s_map(df[[colnames(df)[i]]],
#                           E = E_value, theta = theta_value,
#                           stats_only = FALSE)
#   Smap_results_df <- rbind(Smap_results_df,
#                            data.frame(Variable = colnames(df)[i],
#                                       E_value = smap_res$E,
#                                       theta_value = smap_res$theta,
#                                       rho = simplex_res$rho,
#                                       rmse = simplex_res$rmse))
# }
# write.csv(Smap_results_df, paste0(res_save_path, "parameter.csv"))

## ------------------------------------------------------------------------------ ##
## -------------------------- Seasonal surrogate (UIC) -------------------------- ##
## ------------------------------------------------------------------------------ ##
# ## E test code
# libVar <- "A_H1N1"
# tarVar <- "A_H3N2"
# numSurr <- 2000
# tp_range <- c(-12:4)
# E_range <- c(0:20)
#
# edm_simplex <- rEDM::simplex(df[[libVar]], E = E_range, tp = tp_range, silent = T, stats_only = FALSE)
# edm_BestEforEachTP_res <- list()
# for (tp_var in tp_range) {
#   edm_df_EachTP <- edm_simplex %>% dplyr::filter(tp == tp_var)
#   edm_BestEforEachTP <- edm_df_EachTP[which.min(edm_df_EachTP$rmse),"E"]
#   edm_BestEforEachTP_res[[tp_var + 13]] <- data.frame(edm_df_EachTP %>% dplyr::filter(E == edm_BestEforEachTP))
# }
# edm_BestEforEachTP_res <- do.call(rbind, edm_BestEforEachTP_res)
# # write.csv(BestEforEachTP_res, paste0(uic_surr_tbl_path, "result/simplexE_res_", tarVar, "_cause_", libVar, ".csv"))
# edm_BestE_uic_simplex <- with(edm_BestEforEachTP_res, max(c(0, E)))
# # write.csv(BestE_uic_simplex, paste0(uic_surr_tbl_path, "result/simplex_E_", tarVar, "_cause_", libVar, ".csv"))
# print(paste("EDM - simplex: ", libVar, edm_BestE_uic_simplex))
#
# UIC_simplex_res(libVar, E_range, tp_range)
# UIC_multiSimplex_res(libVar, tarVar, E_range, tp_range)
#
# ## tarVar cause libVar?
# uic_res <- uic.optimal(df, lib_var = libVar, tar_var = tarVar, E = E_range, tau = 1, tp = tp_range, alpha = 0.05)
# print(paste("E - uic optimal: ", libVar, uic_res$E[1] + 1))

# List of influenza subtypes
subtypes <- c("A_H1N1", "A_H3N2", "B")

# Perform UIC analysis for all pairs of subtypes
for (effect_var in subtypes) {
  for (cause_var in subtypes) {
    if (effect_var != cause_var) {
      numSurr <- 2000
      print(paste("Effect:", effect_var, "| Cause:", cause_var))
      Surr_UIC(effect_var, cause_var, numSurr)
    }
  }
}

# Individual calls if needed
Surr_UIC("A", "B", 2000)
Surr_UIC("B", "A", 2000)

## ----------------------------------------------------------------------- ##
## -------------------------- Regularized S-map -------------------------- ##
## ----------------------------------------------------------------------- ##
# Only calculate significant pair in UIC
# cause	    effected
# A_H3N2	A_H1N1
# B	        A_H1N1
# B	        A_H3N2

df <- data.frame(Date = dt0$date,
                 A = dt0$A / mean(dt0$A),
                 A_H1N1 = dt0$A_H1N1 / mean(dt0$A_H1N1),
                 A_H3N2 = dt0$A_H3N2 / mean(dt0$A_H3N2),
                 B = dt0$B / mean(dt0$B)
)
# Perform regularized S-map analysis for specific pairs
# tar_var(cause var.) cause lib_var(effected var)?
# reg_smap_func("B", "A")            # A causes B?

reg_smap_func("A_H1N1", "A_H3N2")  # A/H3N2 causes A/H1N1?
reg_smap_func("A_H1N1", "B")       # B causes A/H1N1?
reg_smap_func("A_H3N2", "B")       # B causes A/H3N2?
