rm(list = ls())
library(remotes)
# install.packages("devtools")
remotes::install_github("ha0ye/rEDM")
remotes::install_github("yutakaos/rUIC")
remotes::install_github("ong8181/macamts")
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
  rUIC,
  sinaplot,
  ggExtra,
  ggdensity
)
# Load library
packageVersion("rEDM") # v 0.7.5
packageVersion("macamts") # v 0.1.4/ v 0.2.0 2025/08/14
packageVersion("rUIC") # v 0.9.12/ v 0.9.15 2025/08/14
# library(rEDM)
# vignette("rEDM-tutorial")
theme_set(theme_cowplot())

## Set path
setwd("/Users/ayo/Desktop/_GSAIS_/mzmtlab/microbiome dynamics/")
## result save path
global_res_save_path <- sprintf("result/FluSub_JP/%s/", Sys.Date())
dir.create(file.path(global_res_save_path))
print(global_res_save_path) # check the path

## s-map coef/pred save path
coef_save_path <- paste0(res_save_path, "smap_coef/")
pred_save_path <- paste0(res_save_path, "smap_predRes/")

## s-map coef/pred save path
coef_save_path <- paste0(res_save_path, "smap_coef/")
pred_save_path <- paste0(res_save_path, "smap_predRes/")
dir.create(file.path(coef_save_path))
dir.create(file.path(pred_save_path))

## plot save path
fig_date_path <- paste0(res_save_path, "figure/")
dir.create(file.path(fig_date_path))
# Raw time series
raw_ts_path <- paste0(fig_date_path, "raw_ts/")
dir.create(file.path(raw_ts_path))
# Standalized time series
standalized_ts_path <- paste0(fig_date_path, "standarlized_incidence_H1H3B.tiff")
dir.create(file.path(standalized_ts_path))
# UIC figure
uic_plot_path <- paste0(res_save_path, "figure/uicFig/")
dir.create(file.path(uic_plot_path))
xmap_plot_path <- paste0(res_save_path, "figure/crossMapping/")
dir.create(file.path(xmap_plot_path))
# tables
dir.create(file.path(paste0(res_save_path, "table/")))
uic_surr_tbl_path <- paste0(res_save_path, "table/uic/")
dir.create(file.path(uic_surr_tbl_path))


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
                 # A = log(dt0$A + 1),
                 B_Victoria = log(dt0$B_Victoria + 1),
                 B_Yamagata = log(dt0$B_Yamagata + 1),

                 A_H1N1 = log(dt0$A_H1N1 + 1),
                 A_H3N2 = log(dt0$A_H3N2 + 1),
                 B = log(dt0$B + 1)
)

## --------------------------------------------------------------- ##
## -------------------------- MDR S-map -------------------------- ##
## --------------------------------------------------------------- ##

# # data_4sp_std <- as.data.frame(apply(data_4sp, 2, function(x) as.numeric(scale(x))))
# # effect_var <- "Trachurus.japonicus"
# effected_var <- "A_H1N1"
# tp_range <- c(-12:0)
# E_range <- c(0:20)
# # Step. 1: Estimate optimal embeding dimension
# simp_x <- rUIC::simplex(df, lib_var = effected_var, E = E_range, tp = 1)
# (Ex <- simp_x[which.min(simp_x$rmse),"E"])
#
# # Step 2: Perform UIC to detect causality
# uic_res <- uic_across(df[,2:ncol(df)], effected_var, E_range = E_range, tp_range = tp_range, silent = TRUE)
#
# # Step 3: Make block to calculate multiview distance
# block_mvd <- make_block_mvd(df[,2:ncol(df)], uic_res, effected_var, E_effect_var = Ex, include_var = "strongest_only", p_threshold = 0.05)
#
# # Step. 4: Compute multiview distance
# multiview_dist <- compute_mvd(block_mvd, effected_var, E = Ex, tp = 1)
# #as.matrix(dist(block_mvd))
#
# # Step. 5: Do MDR S-map
# mdr_res <- s_map_mdr(block_mvd,
#                      dist_w = multiview_dist,
#                      #dist_w = as.matrix(dist(block_mvd)),
#                      theta = 1, # Check!
#                      tp = 1,
#                      regularized = FALSE,
#                      lambda = 0, # Check!
#                      save_smap_coefficients = TRUE)
# mdr_res$stats

## ----------------------------------------------------------------------- ##
## -------------------------- Regularized S-map -------------------------- ##
## ----------------------------------------------------------------------- ##
library("macamts")
# tar_var to lib_var?
effected_var <- "A_H1N1" # effected
cause_var <- c("A_H3N2", "B") # cause

## calculate the Best E
tp_range <- c(-12:4)
E_range <- c(0:20)
## tarVar cause libVar?
simp_effected <- rUIC::simplex(df, lib_var = effected_var, E = E_range, tau=1, tp=1, alpha=0.05)
BestE_effected <- simp_effected[which.min(simp_effected$rmse),"E"]
lag_block_target <- make_block(df[[cause_var[1]]], max_lag = BestE_effected)
print(paste("E - uic simplex: ", lib_var, BestE_effected))


for (var in 1:length(lib_var)) {
  uic_res <- read.csv(paste0(uic_surr_tbl_path, tar_var, "_cause_", lib_var[[var]], "_uic_p_surr_res.csv")) # use uic saving path
  uic_res_pval <- uic_res[uic_res$pval < 0.05,]
  BestTP <- uic_res_pval[which.max(uic_res_pval$ete),]$tp

  if (BestTP == 0) {
    lag_block_lib <- dplyr::lead(df[[lib_var]], n = 1) # tp = 0
  }else {
    lag_block_lib <- dplyr::lag(df[[lib_var]], n = abs(BestTP) - 1)
  }
}


smap_block <- cbind(lag_block_target,
                    lag_block_lib)
smap_block <- smap_block[, 2:ncol(smap_block)]

theta_range <- c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
stat_res <- list()
for (theta in theta_range) {
  print(theta)
  rsmap_ridge <- macamts::extended_lnlp(smap_block, theta = theta, lambda = 0.1,
                                        regularized = T, alpha = 0, random_seed = 1234, save_smap_coefficients = TRUE)
  stat_res <- rbind(stat_res, rsmap_ridge$stats)
}
stat_res <- cbind(theta_range, stat_res)
BestTheta <- stat_res[which.min(stat_res$rmse), "theta_range"]
print(paste("theta - reg. s-map: ", lib_var, tar_var, BestTheta))

rsmap_ridge <- macamts::extended_lnlp(smap_block, theta = BestTheta, lambda = 0.1,
                                      regularized = T, alpha = 0, random_seed = 1234, save_smap_coefficients = TRUE)
smap_pred_res <- rsmap_ridge$model_output
coef_res <- rsmap_ridge$smap_coefficients

maxValue <- max(max(na.omit(smap_pred_res$obs)), max(na.omit(smap_pred_res$pred))) + 1
ggplot(smap_pred_res, aes(x = obs, y = pred, color = time)) +
  geom_abline(slope = 1, linetype = "dashed", color = "black") +
  geom_point() +
  # same xlim, ylim
  xlim(c(NA, maxValue)) +
  ylim(c(NA, maxValue))

coef_res_cols <- data.frame(time = coef_res$time,
                            c_1 = coef_res$c_1,
                            c_2 = coef_res[[paste0("c_", BestE + 1)]])

c1_var <- tar_var
c1_math_expression <- substitute(frac(partialdiff * x1, partialdiff * x2), list(x1 = bquote(.(as.symbol(c1_var))(t + 1)), x2 = bquote(.(as.symbol(c1_var))(t))))
c2_var <- lib_var
c2_math_expression <- substitute(frac(partialdiff * x1, partialdiff * x2), list(x1 = bquote(.(as.symbol(c1_var))(t + 1)), x2 = bquote(.(as.symbol(c2_var))(t))))


coef_res_long <- gather(coef_res, key = "coef", value = "value", c_1, c_2)

# sinaplot, boxplot
ggplot(coef_res_long, aes(x = coef, y = value)) +
  geom_boxplot() +
  geom_violin(alpha = .5) +
  geom_sina(alpha = .5) +

  scale_x_discrete(labels = c(
    "c_1" = c1_math_expression,
    "c_2" = c2_math_expression
  )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  ## mean
  geom_hline(yintercept = mean(na.omit(coef_res$c_1)), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = mean(na.omit(coef_res$c_2)), linetype = "dashed", color = "gray") +

  theme(text = element_text(size = 20), legend.position = c(.25, .65), axis.title.x = element_blank()) +
  labs(x = paste0("Interactions Between ", c1_var, " and ", c2_var, y = "S-map coefficients"))

date_xaxis <- as.POSIXct(df$Date, format = "%Y-%m-%d")

# ts of observed and predicted
ObsPred_plt <- ggplot(smap_pred_res, aes(date_xaxis)) +
  geom_line(aes(y = obs, linetype = "Observation")) +
  geom_line(aes(y = pred, linetype = "Prediction")) +
  scale_colour_grey(start = 0, end = .8) +
  scale_x_datetime(
    date_breaks = "1 year",
    date_labels = "%Y",
    minor_breaks = "1 month"
  ) +
  ylab("Incidence") +
  labs(linetype = c1_var) +
  theme(text = element_text(size = 20), axis.title.x = element_blank(), legend.position = c(.80, .7),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
ObsPred_plt

# plot of c_1
c1Y_plt <- ggplot(coef_res, aes(date_xaxis)) +
  geom_hline(yintercept = quantile(na.omit(coef_res$c_1), 0.25), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = quantile(na.omit(coef_res$c_1), 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = quantile(na.omit(coef_res$c_1), 0.75), linetype = "dashed", color = "gray") +
  geom_line(aes(y = c_1)) +
  # geom_line(aes(y = c_2, linetype = "c_2")) +
  scale_colour_grey(start = 0, end = .8) +
  scale_x_datetime(
    date_breaks = "1 year",
    date_labels = "%Y",
    minor_breaks = "1 month"
  ) +
  labs(y = c1_math_expression) +
  theme(text = element_text(size = 20), axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
c1Y_plt

# plot of c_2
c2X_plt <- ggplot(coef_res, aes(date_xaxis)) +
  geom_hline(yintercept = quantile(na.omit(coef_res$c_2), 0.25), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = quantile(na.omit(coef_res$c_2), 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = quantile(na.omit(coef_res$c_2), 0.75), linetype = "dashed", color = "gray") +
  geom_line(aes(y = c_2)) +
  scale_colour_grey(start = 0, end = .8) +
  labs(y = c2_math_expression) +
  scale_x_datetime(
    date_breaks = "1 year",
    date_labels = "%Y",
    minor_breaks = "1 month"
  ) +
  # geom_hline(yintercept = cor_result, linetype = "dashed", color = "red") +
  theme(text = element_text(size = 20), axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
c2X_plt

# merge 3 plots
merged_plt <- plot_grid(ObsPred_plt, c1Y_plt, c2X_plt, align = "v", cols = 1)
merged_plt
