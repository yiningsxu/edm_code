###########################
## This program is for "Exploring the interaction of influenza subtypes H1N1, H3N2 and B based on empirical dynamics modeling"
## By Yining XU
## latest update: 2024/09/25
## Checking the data/saving path is required before run this code.
## It's necessary to edit the vars while plot time series, bestE, and nonlinear para. theta.(function is being considered)
###########################
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
  glue,
  stats,
  dplyr,
  gridExtra,
  cowplot,
  rlang,
  # macam,
  # rUIC,
  # cowplot,
  ggExtra,
  ggdensity
)
# Load library
packageVersion("rEDM") # v0.7.5
packageVersion("macam") # v0.1.4
packageVersion("rUIC") # v0.9.12
# library(rEDM)
# vignette("rEDM-tutorial")
theme_set(theme_cowplot())

## Set path
setwd("/Users/ayo/Desktop/_GSAIS_/mzmtlab/microbiome dynamics/")
## result save path
res_save_path <- sprintf("result/FluSub_JP/%s/", Sys.Date())

## s-map coef/pred save path
coef_save_path <- paste0(res_save_path,"smap_coef/")
pred_save_path <- paste0(res_save_path,"smap_predRes/")

## plot save path
fig_date_path <- paste0(res_save_path,"figure/")
# Raw time series
raw_ts_path <- paste0(fig_date_path, "raw_ts/")
# Standalized time series
standalized_ts_path <- paste0(fig_date_path, "standarlized_incidence_H1H3B.tiff")
# Best E
bestE_fig_path <- paste0(fig_date_path, "BestE/")
# nonlinear theta
theta_fig_path <- paste0(fig_date_path, "theta/")
# attractor
attractor_path <- paste0(fig_date_path, "FluSub_attractor/")
# interaction boxplot
boxplot_path <- paste0(fig_date_path, "FluSub_interaction/boxplot/")
# interaction time series
ts_path <- paste0(fig_date_path, "FluSub_interaction/ts/")
# ccm rho
ccm_rho_fig_path <- paste0(fig_date_path, "FluSub_causality/rho_")
# ccm rmse
ccm_rmse_fig_path <- paste0(fig_date_path, "FluSub_causality/rmse_")
## --------------------- Preparation --------------------- ##
## import data
# dt0 <- read.csv("data/FluSub_jp/FluSub_10to19_jp.csv")
# dt0 <- dt0[c(70:537),]

dt0 <- read.csv("data/FluSub_jp/FluSub_11to19_jp_per_20240925.csv")

# dt <- subset(dt, select = -c(X))
dt0 <- dt0 %>%
  mutate(date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1")))
dt0$date <- as.Date(dt0$date)
# dt0 <- dt0[, c(9, 3:8)]
# 416 x 14
# dt0$A <- dt0$AH1pdm + dt0$AH3

# dt0[dt0$year == 2009,]
# for (i in 0:9) {
#   print(2010 + i)
#   print(i * 52 + 1)
#   print(i * 52 + 52)
#   start <- i * 52 + 1
#   end <- i * 52 + 52
#   dt_epidYear <- dt0[start:end,]
#   year <- 2009+i
#   dt_epidYear <- dt0[dt0$year == year,]
#   ggplot(dt_epidYear, aes(week)) +
#     # geom_point(aes(y = AH1pdm, colour = "A(H1N1)")) +
#     # geom_point(aes(y = AH3, colour = "A(H3N2)")) +
#     # geom_point(aes(y = B, colour = "B")) +
#
#     geom_line(aes(y = AH1pdm, colour = "A(H1N1)"),linewidth = 1, show.legend = FALSE) +
#     geom_line(aes(y = AH3, colour = "A(H3N2)"),linewidth = 1, show.legend = FALSE) +
#     geom_line(aes(y = B, colour = "B"), linewidth = 1, show.legend = FALSE) +
#
#     labs(colour = "Influenza Subtypes") +
#     scale_colour_grey(start = 0, end = .8) +
#     scale_x_continuous(breaks = seq(0, 52, by = 10)) +
#     ylim(c(0, 1690))+
#     theme_bw() +
#     theme(title = element_text(paste(2009 + i, "-", 2009 + i + 1), size = 20),
#           text = element_text(size = 20),
#           axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none",
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")
#     ) +
#     xlab("Week") +
#     ylab("Incidence")
#   ggsave(paste0("result/figure/incidence_ts/incidence_", 2009 + i, ".tiff"), units = "in", width = 3, height = 3, dpi = 300, compression = 'lzw')
# }
# dt_epidYear <- dt0[521:537,]
# ggplot(dt0, aes(week)) +
#   # geom_point(aes(y = AH1pdm, colour = "A(H1N1)")) +
#   # geom_point(aes(y = AH3, colour = "A(H3N2)")) +
#   # geom_point(aes(y = B, colour = "B")) +
#
#   geom_line(aes(y = AH1pdm, colour = "A(H1N1)"),linewidth = 1) +
#   geom_line(aes(y = AH3, colour = "A(H3N2)"),linewidth = 1) +
#   geom_line(aes(y = B, colour = "B"), linewidth = 1) +
#
#   labs(colour = "Influenza Subtypes") +
#   scale_colour_grey(start = 0, end = .8) +
#   # scale_x_continuous(breaks = seq(0, nrow(dt_epidYear), by = 10)) +
#   theme_bw() +
#   # theme(title = element_text(paste(2009 + i, "-", 2009 + i + 1), size = 20)) +
#   # theme(text = element_text(size = 20),
#   #       axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none",
#   #       panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#   facet_grid(year~.) +
#   xlab("Epidemiological Week") +
#   ylab("Incidence")
# ggsave(paste0("result/figure/incidence_ts/incidence_2019_w36-52.tiff"), units = "in", width = 8, height = 8, dpi = 300, compression = 'lzw')

## ----------------- Figure Raw time series ----------------- ##
ggplot(dt0, aes(date)) +
  geom_point(aes(y = A, colour = "A")) +
  geom_point(aes(y = B, colour = "B")) +

  geom_line(aes(y = A, colour = "A")) +
  geom_line(aes(y = B, colour = "B")) +

  labs(colour = "Influenza Subtypes") +
  scale_colour_grey(start = 0, end = .6) +
  theme_bw() +
  ylim(c(0, 1600)) +
  theme(text = element_text(size = 20)) +
  xlab("Date") +
  ylab("Incidence")
ggsave(paste0(raw_ts_path, "incidence_AB.tiff"), units = "in", width = 16, height = 6, dpi = 300, compression = 'lzw')

ggplot(dt0, aes(date)) +
  geom_point(aes(y = AH1pdm, colour = "A/H1N1")) +
  geom_point(aes(y = AH3, colour = "A/H3N2")) +

  geom_line(aes(y = AH1pdm, colour = "A/H1N1")) +
  geom_line(aes(y = AH3, colour = "A/H3N2")) +

  labs(colour = "Influenza Subtypes") +
  scale_colour_grey(start = 0, end = .6) +
  theme_bw() +
  ylim(c(0, 1600)) +
  theme(text = element_text(size = 20)) +
  xlab("Date") +
  ylab("Incidence")
ggsave(paste0(raw_ts_path, "incidence_AH1H3.tiff"), units = "in", width = 16, height = 6, dpi = 300, compression = 'lzw')

ggplot(dt0, aes(date)) +
  # geom_point(aes(y = B, colour = "B(Total)")) +
  geom_point(aes(y = Bvic, colour = "B/Victoria")) +
  geom_point(aes(y = Byama, colour = "B/Yamagata")) +

  # geom_line(aes(y = B, colour = "B(Total)")) +
  geom_line(aes(y = Bvic, colour = "B/Victoria")) +
  geom_line(aes(y = Byama, colour = "B/Yamagata")) +

  labs(colour = "Influenza Subtypes") +
  scale_colour_grey(start = 0, end = .6) +
  theme_bw() +
  ylim(c(0, 1600)) +
  theme(text = element_text(size = 20)) +
  xlab("Date") +
  ylab("Incidence")
ggsave(paste0(raw_ts_path, "incidence_BvicByama.tiff"), units = "in", width = 16, height = 6, dpi = 300, compression = 'lzw')


#
dt <- dt0[, c(1, 8, 4, 2:3, 5:7)]
dt_stad <- data.frame(scale(dt[, 2:length(dt)]))
dt_stad <- mutate(dt_stad, Date = dt$date)
dt_stad <- dt_stad[, c(length(dt_stad), 1:length(dt_stad) - 1)]
write.csv(dt_stad, "data/FluSub_jp/standarlized_FluSub_10to19_jp.csv")
#
# dt_norm <- dt
# for (i in colnames(dt_norm)[2:length(dt_norm)]){
#   dt_norm[[i]] <- dt_norm[[i]]/mean(dt_norm[[i]])
# }
# write.csv(dt_norm,"data/FluSub_jp/normalized_FluSub_10to19_jp.csv")

## --------------------- Analysis --------------------- ##
## import data
# dt <- read.csv("data/FluSub_jp/normalized_FluSub_10to19_jp.csv")
dt <- read.csv("data/FluSub_jp/standarlized_FluSub_10to19_jp.csv")
dt <- subset(dt, select = -c(X))
# 468 x 7
dt$Date <- as.Date(dt$Date)
str(dt)
# dt <- dt[, c(1:4 )]
colnames(dt)
# [1] "Date"   "AH1pdm" "AH3"
# [4] "B"      "Bvic"   "Byama"
# [7] "Bunk"

## ----------------- Figure Standarlized time series ----------------- ##
dt_plot <- ggplot(dt, aes(Date)) +
  geom_point(aes(y = AH1pdm, colour = "A(H1N1)")) +
  geom_point(aes(y = AH3, colour = "A(H3N2)")) +
  geom_point(aes(y = B, colour = "B")) +
  # geom_point(aes(y = Bvic, colour = "Bvic")) +
  # geom_point(aes(y = Byama, colour = "Byama")) +
  # geom_point(aes(y = Bunk, colour = "Bunk")) +

  geom_line(aes(y = AH1pdm, colour = "A(H1N1)")) +
  geom_line(aes(y = AH3, colour = "A(H3N2)")) +
  geom_line(aes(y = B, colour = "B")) +
  # geom_line(aes(y = Bvic, colour = "Bvic")) +
  # geom_line(aes(y = Byama, colour = "Byama")) +
  # geom_line(aes(y = Bunk, colour = "Bunk")) +

  labs(colour = "Influenza Subtypes") +
  scale_colour_grey(start = 0, end = .8) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  xlab("Date") +
  ylab("Standardized Incidence")
dt_plot

ggsave(paste0(standalized_ts_path), units = "in", width = 16, height = 8, dpi = 300, compression = 'lzw')

## --------------------- Simple Projection --------------------- ##
# 最適埋め込み次元の推定
simplex_res_bestE <- data.frame(matrix(nrow = 10))
# simplex_res_bestE[["E"]] <- rEDM::simplex(dt$AH1pdm, E = 1:20,silent = T)$E
simplex_res_bestE[["E"]] <- rEDM::simplex(dt$AH1pdm, E = 1:10, silent = T)$E
for (i in 2:length(dt)) {
  simplex_res <- rEDM::simplex(dt[[colnames(dt)[i]]],
                               E = 1:10,
                               silent = T,
                               stats_only = FALSE)$rmse
  simplex_res_bestE[[colnames(dt)[i]]] <- simplex_res
}
sum_simp_res <- simplex_res_bestE[, c(2:length(simplex_res_bestE))]

# res_embedE <- EmbedDimension(dataFrame = dt, lib = "1 260", pred = "1 260", columns = "AH1pdm", target = "AH1pdm", showPlot = TRUE)

ggplot(sum_simp_res, aes(x = E)) +
  geom_point(aes(y = A)) +
  geom_point(aes(y = B)) +
  geom_point(aes(y = AH1pdm)) +
  geom_point(aes(y = AH3)) +
  geom_point(aes(y = Bvic)) +
  geom_point(aes(y = Byama)) +
  # geom_point(aes(y = Bunk)) +

  geom_line(aes(y = A, linetype = "A")) +
  geom_line(aes(y = B, linetype = "B")) +
  geom_line(aes(y = AH1pdm, linetype = "A/H1N1")) +
  geom_line(aes(y = AH3, linetype = "A/H3N2")) +
  geom_line(aes(y = Bvic, linetype = "B/Victoria")) +
  geom_line(aes(y = Byama, linetype = "B/Yamagata")) +
  # geom_line(aes(y = Bunk, linetype = "Bunk")) +

  geom_point(aes(x = sum_simp_res[which.min(sum_simp_res$A),]$E,
                 y = sum_simp_res[which.min(sum_simp_res$A),]$A),
             shape = 18, size = 5) +
  geom_point(aes(x = sum_simp_res[which.min(sum_simp_res$B),]$E,
                 y = sum_simp_res[which.min(sum_simp_res$B),]$B),
             shape = 18, size = 5) +
  geom_point(aes(x = sum_simp_res[which.min(sum_simp_res$AH1pdm),]$E,
                 y = sum_simp_res[which.min(sum_simp_res$AH1pdm),]$AH1pdm),
             shape = 18, size = 5) +
  geom_point(aes(x = sum_simp_res[which.min(sum_simp_res$AH3),]$E,
                 y = sum_simp_res[which.min(sum_simp_res$AH3),]$AH3),
             shape = 18, size = 5) +
  geom_point(aes(x = sum_simp_res[which.min(sum_simp_res$Bvic),]$E,
                 y = sum_simp_res[which.min(sum_simp_res$Bvic),]$Bvic),
             shape = 18, size = 5) +
  geom_point(aes(x = sum_simp_res[which.min(sum_simp_res$Byama),]$E,
                 y = sum_simp_res[which.min(sum_simp_res$Byama),]$Byama),
             shape = 18, size = 5) +
  # geom_point(aes(x = sum_simp_res[which.min(sum_simp_res$Bunk),]$E,
  #                y = sum_simp_res[which.min(sum_simp_res$Bunk),]$Bunk),
  #            shape = 18, size = 5) +

  scale_x_continuous(breaks = seq(1, 20, by = 1)) +
  labs(linetype = "Influenza Subtypes") +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  # , legend.position = c(.65, .85)
  xlab("Embedding Dimension(E)") +
  ylab("Forecaset Skill (RMSE)")

ggsave(paste0(bestE_fig_path, "BestE_all.tiff"), units = "in", width = 10, height = 8, dpi = 300, compression = 'lzw')

## --------------------- Result of Simplex Projection --------------------- ##
simplex_res_list <- list()
for (i in 2:length(dt)) {
  E_value <- sum_simp_res[which.min(sum_simp_res[[colnames(dt)[i]]]), "E"]
  simplex_res <- rEDM::simplex(dt[[colnames(dt)[i]]],
                               E = E_value,
                               stats_only = FALSE)
  simplex_res_list[[colnames(dt)[i]]] <- list(simplex_res)
  print(colnames(dt)[i])
  print(paste("rho --------------------", simplex_res$E))
  print(paste("rho --------------------", simplex_res$rho))
  print(paste("rmse -------------------", simplex_res$rmse))
}
simplex_res_list

### --------------------------- S-map --------------------------- ##
uniSmap_res_list <- list()
for (i in 2:length(dt)) {
  E_value <- sum_simp_res[which.min(sum_simp_res[[colnames(dt)[i]]]), "E"]
  uniSmap_res <- s_map(dt[[colnames(dt)[i]]],
                       E = E_value,
                       silent = T,
                       stats_only = FALSE)
  uniSmap_res_list[[colnames(dt)[i]]] <- list(uniSmap_res)
}

sum_smap_res <- data.frame("theta" = data.frame(uniSmap_res_list[["AH1pdm"]])$theta)
for (i in 2:ncol(dt)) {
  print(colnames(dt)[i])
  rmse_value <- data.frame(uniSmap_res_list[[colnames(dt)[i]]])$rmse
  rmse <- data.frame(rmse_value)
  colnames(rmse) <- colnames(dt)[i]
  sum_smap_res <- cbind(sum_smap_res, rmse)
}

smap_res_list <- list()
for (i in 2:length(dt)) {
  E_value <- sum_simp_res[which.min(sum_simp_res[[colnames(dt)[i]]]), "E"]
  theta_value <- sum_smap_res[which.min(sum_smap_res[[colnames(dt)[i]]]), "theta"]
  # theta_value <- 0
  smap_res <- rEDM::s_map(dt[[colnames(dt)[i]]],
                          E = E_value, theta = theta_value,
                          stats_only = FALSE)
  smap_res_list[[colnames(dt)[i]]] <- list(smap_res)
  print(colnames(dt)[i])
  print(paste("E --------------------", smap_res$E))
  print(paste("theta --------------------", smap_res$theta))
  print(paste("rho --------------------", smap_res$rho))
  print(paste("rmse -------------------", smap_res$rmse))
}
smap_res_list

ggplot(sum_smap_res, aes(theta)) +
  geom_point(aes(y = A)) +
  geom_point(aes(y = B)) +
  geom_point(aes(y = AH1pdm)) +
  geom_point(aes(y = AH3)) +
  geom_point(aes(y = Bvic)) +
  geom_point(aes(y = Byama)) +
  # geom_point(aes(y = Bunk)) +

  geom_line(aes(y = A, linetype = "A")) +
  geom_line(aes(y = B, linetype = "B")) +
  geom_line(aes(y = AH1pdm, linetype = "A/H1N1")) +
  geom_line(aes(y = AH3, linetype = "A/H3N2")) +
  geom_line(aes(y = Bvic, linetype = "B/Victoria")) +
  geom_line(aes(y = Byama, linetype = "B/Yamagata")) +

  geom_point(aes(x = sum_smap_res[which.min(sum_smap_res$A),]$theta,
                 y = sum_smap_res[which.min(sum_smap_res$A),]$A),
             shape = 18, size = 5) +
  geom_point(aes(x = sum_smap_res[which.min(sum_smap_res$B),]$theta,
                 y = sum_smap_res[which.min(sum_smap_res$B),]$B),
             shape = 18, size = 5) +
  geom_point(aes(x = sum_smap_res[which.min(sum_smap_res$AH1pdm),]$theta,
                 y = sum_smap_res[which.min(sum_smap_res$AH1pdm),]$AH1pdm),
             shape = 18, size = 5) +
  geom_point(aes(x = sum_smap_res[which.min(sum_smap_res$AH3),]$theta,
                 y = sum_smap_res[which.min(sum_smap_res$AH3),]$AH3),
             shape = 18, size = 5) +

  geom_point(aes(x = sum_smap_res[which.min(sum_smap_res$Bvic),]$theta,
                 y = sum_smap_res[which.min(sum_smap_res$Bvic),]$Bvic),
             shape = 18, size = 5) +
  geom_point(aes(x = sum_smap_res[which.min(sum_smap_res$Byama),]$theta,
                 y = sum_smap_res[which.min(sum_smap_res$Byama),]$Byama),
             shape = 18, size = 5) +
  # geom_point(aes(x = sum_smap_res[which.min(sum_smap_res$Bunk),]$theta,
  #                y = sum_smap_res[which.min(sum_smap_res$Bunk),]$Bunk),
  #            shape = 18, size = 5) +

  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  labs(linetype = "Influenza Subtypes") +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  # theme(legend.position="none")+
  # , legend.position = c(.85, .65)
  xlab("Nonlinear Parameter (θ)") +
  ylab("Forecaset Skill (RMSE)")

ggsave(paste0(theta_fig_path, "theta_all.tiff"), units = "in", width = 10, height = 8, dpi = 300, compression = 'lzw')

### --------------------------- Multivariate S-map --------------------------- ##

multiSmap_func <- function(Smap_colname) {
  # data input to Multivar S-map
  smap_block <- dt[, Smap_colname]
  # print the column name and the name of target column(always column 1)
  print(colnames(smap_block))
  print(paste0("target column : ", colnames(smap_block)[1]))

  # value of thetas
  theta <- c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
  # find the theta (by rmse)
  multiSmap_res <- block_lnlp(smap_block, method = "s-map",
                              theta = theta,
                              target_column = 1,
                              first_column_time = FALSE,
                              silent = TRUE, stats_only = FALSE)
  multiSmap_min_theta <- multiSmap_res[which.min(multiSmap_res$rmse), "theta"]

  # execute Multivariate S-map with the theta
  multiSmap_result <- block_lnlp(smap_block, method = "s-map",
                                 theta = multiSmap_min_theta,
                                 target_column = 1,
                                 first_column_time = FALSE,
                                 silent = TRUE, save_smap_coefficients = TRUE)

  # print the theta, and result of rmse & rho with the theta
  print(paste0("min theta : ", multiSmap_min_theta))
  print(paste0("rmse : ", multiSmap_result$rmse))
  print(paste0("rho : ", multiSmap_result$rho))
  # result of coefficient and prediction result
  coef_res <- data.frame(multiSmap_result$smap_coefficients[[1]])
  smap_pred_res <- data.frame(multiSmap_result$model_output[[1]])
  # write.csv(coef_res, file = paste0(coef_save_path,target_column,".csv"))
  # write.csv(smap_pred_res, file = paste0(pred_save_path,"BAH1AH3.csv"))


  ### plot the result of prediction
  maxValue <- max(max(na.omit(smap_pred_res$obs)), max(na.omit(smap_pred_res$pred))) + 1
  ggplot(smap_pred_res, aes(x = obs, y = pred, color = time)) +
    geom_abline(slope = 1, linetype = "dashed", color = "black") +
    geom_point() +
    # same xlim, ylim
    xlim(c(NA, maxValue)) +
    ylim(c(NA, maxValue))
  # colour faded by time (smaller the time is, earlier in the time series)
  # line of 1:1 ratio
  if (ncol(smap_block) == 2) {
    # save
    ggsave(paste0(fig_date_path, "Smap/prediction/", colnames(smap_block)[1], "(target) & ", colnames(smap_block)[2], ".tiff"),
           units = "in",
           width = 5, height = 5, dpi = 300, compression = 'lzw')
    write.csv(coef_res, file = paste0(coef_save_path, colnames(smap_block)[1], "(target) & ", colnames(smap_block)[2], ".csv"))
    write.csv(smap_pred_res, file = paste0(pred_save_path, colnames(smap_block)[1], "(target) & ", colnames(smap_block)[2], ".csv"))

    # statistic result of coef result
    print("==== mean / SD ====")
    print(paste("coef c1",colnames(smap_block)[1],"->",mean(na.omit(coef_res$c_1)),sd(na.omit(coef_res$c_1))))
    print(paste("coef c2",colnames(smap_block)[2],"->",mean(na.omit(coef_res$c_2)),sd(na.omit(coef_res$c_2))))

  }else {
    ggsave(paste0(fig_date_path, "Smap/prediction/", colnames(smap_block)[1], "(target) & ", colnames(smap_block)[2], "_", colnames(smap_block)[3], ".tiff"),
           units = "in",
           width = 5, height = 5, dpi = 300, compression = 'lzw')
    write.csv(coef_res, file = paste0(coef_save_path, colnames(smap_block)[1], "(target) & ", colnames(smap_block)[2], "_", colnames(smap_block)[3], ".csv"))
    write.csv(smap_pred_res, file = paste0(pred_save_path, colnames(smap_block)[1], "(target) & ", colnames(smap_block)[2], "_", colnames(smap_block)[3], ".csv"))
    # statistic result of coef result
    print("==== mean / SD ====")
    print(paste("coef c1",colnames(smap_block)[1],"->",mean(na.omit(coef_res$c_1)),sd(na.omit(coef_res$c_1))))
    print(paste("coef c2",colnames(smap_block)[2],"->",mean(na.omit(coef_res$c_2)),sd(na.omit(coef_res$c_2))))
    print(paste("coef c3",colnames(smap_block)[3],"->",mean(na.omit(coef_res$c_3)),sd(na.omit(coef_res$c_3))))
  }
  if (ncol(smap_block) == 2) {
    coef_res_long <- gather(coef_res, key = "coef", value = "value", c_1, c_2)

    c1_var <- colnames(smap_block)[1]
    c1_math_expression <- substitute(frac(partialdiff * x1, partialdiff * x2), list(x1 = bquote(.(as.symbol(c1_var))(t + 1)), x2 = bquote(.(as.symbol(c1_var))(t))))
    c2_var <- colnames(smap_block)[2]
    c2_math_expression <- substitute(frac(partialdiff * x1, partialdiff * x2), list(x1 = bquote(.(as.symbol(c1_var))(t + 1)), x2 = bquote(.(as.symbol(c2_var))(t))))

    boxplot <- ggplot(coef_res_long, aes(x = coef, y = value)) +
      geom_boxplot() +
      # geom_point()+
      # geom_hdr_lines(xlim = c(100, 300), ylim = c(0, 100))+
      theme_bw() +
      scale_x_discrete(labels = c(
        "c_1" = c1_math_expression,
        "c_2" = c2_math_expression
      )) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      theme(text = element_text(size = 20), legend.position = c(.25, .65), axis.title.x = element_blank()) +
      labs(x = paste0("Interactions Between ", colnames(smap_block)[1], " and ", colnames(smap_block)[2]), y = "S-map coefficients")
    ggsave(paste0(fig_date_path, "Smap/boxplot/", colnames(smap_block)[1], "(target) & ", colnames(smap_block)[2], ".tiff"),
           plot =boxplot,
           units = "in",
           width = 7, height = 5, dpi = 300, compression = 'lzw')

    date_xaxis <- as.POSIXct(dt$Date, format = "%Y-%m-%d")

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
      labs(linetype = colnames(smap_block)[1]) +
      theme(text = element_text(size = 20), axis.title.x = element_blank(), legend.position = c(.80, .7),
            axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))

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

    merged_plt <- plot_grid(ObsPred_plt, c1Y_plt, c2X_plt, align = "v", cols = 1)
    ggsave(paste0(fig_date_path, "Smap/coef/", colnames(smap_block)[1], "(target) & ", colnames(smap_block)[2], ".tiff"),
           plot = merged_plt,
           units = "in",
           width = 12, height = 9, dpi = 300, compression = 'lzw')

  }else {
    coef_res_long <- gather(coef_res, key = "coef", value = "value", c_1, c_2, c_3)
    c1_var <- colnames(smap_block)[1]
    c1_math_expression <- substitute(frac(partialdiff * x1, partialdiff * x2), list(x1 = bquote(.(as.symbol(c1_var))(t + 1)), x2 = bquote(.(as.symbol(c1_var))(t))))
    c2_var <- colnames(smap_block)[2]
    c2_math_expression <- substitute(frac(partialdiff * x1, partialdiff * x2), list(x1 = bquote(.(as.symbol(c1_var))(t + 1)), x2 = bquote(.(as.symbol(c2_var))(t))))
    c3_var <- colnames(smap_block)[3]
    c3_math_expression <- substitute(frac(partialdiff * x1, partialdiff * x2), list(x1 = bquote(.(as.symbol(c1_var))(t + 1)), x2 = bquote(.(as.symbol(c3_var))(t))))

    boxplot <- ggplot(coef_res_long, aes(x = coef, y = value)) +
      geom_boxplot() +
      # geom_point()+
      # geom_hdr_lines(xlim = c(100, 300), ylim = c(0, 100))+
      theme_bw() +
      scale_x_discrete(labels = c(
        "c_1" = c1_math_expression,
        "c_2" = c2_math_expression,
        "c_3" = c3_math_expression
      )) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      theme(text = element_text(size = 20), legend.position = c(.25, .65), axis.title.x = element_blank(),) +
      labs(x = paste0("Interactions Between ", colnames(smap_block)[1], " and ", colnames(smap_block)[2]), y = "S-map coefficients")
    ggsave(paste0(fig_date_path, "Smap/boxplot/", colnames(smap_block)[1], "(target) & ", colnames(smap_block)[2], "_", colnames(smap_block)[3], ".tiff"),
           plot =boxplot,
           units = "in",
           width = 7, height = 5, dpi = 300, compression = 'lzw')

    date_xaxis <- as.POSIXct(dt$Date, format = "%Y-%m-%d")

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
      labs(linetype = colnames(smap_block)[1]) +
      theme(text = element_text(size = 20), axis.title.x = element_blank(), legend.position = c(.80, .7),
            axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))

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

    c2X1_plt <- ggplot(coef_res, aes(date_xaxis)) +
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
      theme(text = element_text(size = 20), axis.title.x = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))

    c2X2_plt <- ggplot(coef_res, aes(date_xaxis)) +
      geom_hline(yintercept = quantile(na.omit(coef_res$c_3), 0.25), linetype = "dashed", color = "gray") +
      geom_hline(yintercept = quantile(na.omit(coef_res$c_3), 0.5), linetype = "dashed", color = "black") +
      geom_hline(yintercept = quantile(na.omit(coef_res$c_3), 0.75), linetype = "dashed", color = "gray") +
      geom_line(aes(y = c_3)) +
      scale_colour_grey(start = 0, end = .8) +
      labs(y = c3_math_expression) +
      scale_x_datetime(
        date_breaks = "1 year",
        date_labels = "%Y",
        minor_breaks = "1 month"
      ) +
      theme(text = element_text(size = 20), axis.title.x = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))

    merged_plt <- plot_grid(ObsPred_plt, c1Y_plt, c2X1_plt, c2X2_plt, align = "v", cols = 1)
    # save
    ggsave(paste0(fig_date_path, "Smap/coef/", colnames(smap_block)[1], "(target) & ", colnames(smap_block)[2], "_", colnames(smap_block)[3], ".tiff"),
           plot = merged_plt,
           units = "in",
           width = 12, height = 12, dpi = 300, compression = 'lzw')

  }

}

Smap_colname <- c("A", "B")
multiSmap_func(Smap_colname)
Smap_colname <- c("B", "A")
multiSmap_func(Smap_colname)
Smap_colname <- c("AH1pdm", "AH3")
multiSmap_func(Smap_colname)
Smap_colname <- c("AH3", "AH1pdm")
multiSmap_func(Smap_colname)
Smap_colname <- c("B", "AH3")
multiSmap_func(Smap_colname)
Smap_colname <- c("AH3", "B")
multiSmap_func(Smap_colname)
Smap_colname <- c("B", "AH1pdm")
multiSmap_func(Smap_colname)
Smap_colname <- c("AH1pdm", "B")
multiSmap_func(Smap_colname)
Smap_colname <- c("Byama", "Bvic")
multiSmap_func(Smap_colname)
Smap_colname <- c("Bvic", "Byama")
multiSmap_func(Smap_colname)

Smap_colname <- c("AH1pdm", "AH3", "B")
multiSmap_func(Smap_colname)
Smap_colname <- c("AH3", "AH1pdm", "B")
multiSmap_func(Smap_colname)
Smap_colname <- c("B", "AH1pdm", "AH3")
multiSmap_func(Smap_colname)

## --------------------------- CCM --------------------------- ##
ccm_stats_func <- function(ccm_res) {
  # mean
  means <- ccm_means(ccm_res)
  # rho Q1, Q3
  quantile_res_rho <- aggregate(rho ~ lib_size, data = ccm_res, FUN = "quantile", probs = c(25, 75) / 100)
  ccm_quantile_res_rho <- data.frame(quantile_res_rho[[2]])
  colnames(ccm_quantile_res_rho) <- c("rhoQ1", "rhoQ3")
  # rho Q1, Q3
  quantile_res_rmse <- aggregate(rmse ~ lib_size, data = ccm_res, FUN = "quantile", probs = c(25, 75) / 100)
  ccm_quantile_res_rmse <- data.frame(quantile_res_rmse[[2]])
  colnames(ccm_quantile_res_rmse) <- c("rmseQ1", "rmseQ3")

  quantile_res <- cbind(ccm_quantile_res_rho, ccm_quantile_res_rmse)

  ccm_stats <- cbind(means, quantile_res)

  return(ccm_stats)
}

ccm_func <- function(ccm_colname) {
  ccm_block <- dt[, ccm_colname]
  x <- colnames(ccm_block[1])
  y <- colnames(ccm_block[2])

  Ex <- sum_simp_res[which.min(sum_simp_res[[x]]),]$E
  Ey <- sum_simp_res[which.min(sum_simp_res[[y]]),]$E

  ccm_res_xtoy_Ex <- ccm(ccm_block, E = Ex, lib = c(1, NROW(ccm_block)), lib_column = 1, target_column = 2, lib_sizes = seq(Ex + 1, nrow(dt), by = 20))
  ccm_stats_xtoy_Ex <- ccm_stats_func(ccm_res_xtoy_Ex)

  ccm_res_xtoy_Ey <- ccm(ccm_block, E = Ey, lib = c(1, NROW(ccm_block)), lib_column = 1, target_column = 2, lib_sizes = seq(Ey + 1, nrow(dt), by = 20))
  ccm_stats_xtoy_Ey <- ccm_stats_func(ccm_res_xtoy_Ey)

  ccm_res_ytox_Ex <- ccm(ccm_block, E = Ex, lib = c(1, NROW(ccm_block)), lib_column = 2, target_column = 1, lib_sizes = seq(Ex + 1, nrow(dt), by = 20))
  ccm_stats_ytox_Ex <- ccm_stats_func(ccm_res_ytox_Ex)

  ccm_res_ytox_Ey <- ccm(ccm_block, E = Ey, lib = c(1, NROW(ccm_block)), lib_column = 2, target_column = 1, lib_sizes = seq(Ey + 1, nrow(dt), by = 20))
  ccm_stats_ytox_Ey <- ccm_stats_func(ccm_res_ytox_Ey)

  sum_ccm_res_Ex <- data.frame(
    libsize = ccm_stats_ytox_Ex$lib_size,

    YpredictX_rho = ccm_stats_ytox_Ex$rho,
    YpredictX_rmse = ccm_stats_ytox_Ex$rmse,
    YpredictX_Q1_rho = ccm_stats_ytox_Ex$rhoQ1,
    YpredictX_Q3_rho = ccm_stats_ytox_Ex$rhoQ3,
    YpredictX_Q1_rmse = ccm_stats_ytox_Ex$rmseQ1,
    YpredictX_Q3_rmse = ccm_stats_ytox_Ex$rmseQ3,

    XpredictY_rho = ccm_stats_xtoy_Ex$rho,
    XpredictY_rmse = ccm_stats_xtoy_Ex$rmse,
    XpredictY_Q1_rho = ccm_stats_xtoy_Ex$rhoQ1,
    XpredictY_Q3_rho = ccm_stats_xtoy_Ex$rhoQ3,
    XpredictY_Q1_rmse = ccm_stats_xtoy_Ex$rmseQ1,
    XpredictY_Q3_rmse = ccm_stats_xtoy_Ex$rmseQ3
  )

  sum_ccm_res_Ey <- data.frame(
    libsize = ccm_stats_ytox_Ey$lib_size,

    YpredictX_rho = ccm_stats_ytox_Ey$rho,
    YpredictX_rmse = ccm_stats_ytox_Ey$rmse,
    YpredictX_Q1_rho = ccm_stats_ytox_Ey$rhoQ1,
    YpredictX_Q3_rho = ccm_stats_ytox_Ey$rhoQ3,
    YpredictX_Q1_rmse = ccm_stats_ytox_Ey$rmseQ1,
    YpredictX_Q3_rmse = ccm_stats_ytox_Ey$rmseQ3,

    XpredictY_rho = ccm_stats_xtoy_Ey$rho,
    XpredictY_rmse = ccm_stats_xtoy_Ey$rmse,
    XpredictY_Q1_rho = ccm_stats_xtoy_Ey$rhoQ1,
    XpredictY_Q3_rho = ccm_stats_xtoy_Ey$rhoQ3,
    XpredictY_Q1_rmse = ccm_stats_xtoy_Ey$rmseQ1,
    XpredictY_Q3_rmse = ccm_stats_xtoy_Ey$rmseQ3
  )

  ggplot(sum_ccm_res_Ex, aes(libsize)) +
    geom_ribbon(aes(ymin = YpredictX_Q1_rho, ymax = YpredictX_Q3_rho), fill = "grey80", alpha = 0.5) +
    geom_ribbon(aes(ymin = XpredictY_Q1_rho, ymax = XpredictY_Q3_rho), fill = "grey80", alpha = 0.5) +
    geom_line(aes(y = YpredictX_rho, linetype = paste(y, "xmap", x))) +
    geom_line(aes(y = XpredictY_rho, linetype = paste(x, "xmap", y))) +
    geom_point(aes(y = YpredictX_rho)) +
    geom_point(aes(y = XpredictY_rho)) +

    labs(linetype = "Cross Map", subtitle = paste("Optimal Embedding Dimension of ", x, ":", Ex)) +
    theme_bw() +
    scale_x_continuous(breaks = seq(0, 600, by = 100)) +
    theme(text = element_text(size = 20), legend.position = c(.65, .65), plot.subtitle = element_text(size = 12)) +
    xlab("Library Size") +
    ylab("Forecaset Skill (ρ)")
  ggsave(paste0(ccm_rho_fig_path, x, " & ", y, "_BestE(", x, ")_", Ex, ".tiff"),
         units = "in",
         width = 8, height = 8, dpi = 300, compression = 'lzw')

  ggplot(sum_ccm_res_Ex, aes(libsize)) +
    geom_ribbon(aes(ymin = YpredictX_Q1_rmse, ymax = YpredictX_Q3_rmse), fill = "grey80", alpha = 0.5) +
    geom_ribbon(aes(ymin = XpredictY_Q1_rmse, ymax = XpredictY_Q3_rmse), fill = "grey80", alpha = 0.5) +
    geom_line(aes(y = YpredictX_rmse, linetype = paste(y, "xmap", x))) +
    geom_line(aes(y = XpredictY_rmse, linetype = paste(x, "xmap", y))) +
    geom_point(aes(y = YpredictX_rmse)) +
    geom_point(aes(y = XpredictY_rmse)) +

    labs(linetype = "Cross Map", subtitle = paste("Optimal Embedding Dimension of ", x, ":", Ex)) +
    theme_bw() +
    scale_x_continuous(breaks = seq(0, 600, by = 100)) +
    theme(text = element_text(size = 20), legend.position = c(.65, .65), plot.subtitle = element_text(size = 12)) +
    xlab("Library Size") +
    ylab("Forecaset Skill (RMSE)")
  ggsave(paste0(ccm_rmse_fig_path, x, " & ", y, "_BestE(", x, ")_", Ex, ".tiff"),
         units = "in",
         width = 8, height = 8, dpi = 300, compression = 'lzw')

  ggplot(sum_ccm_res_Ey, aes(libsize)) +
    geom_ribbon(aes(ymin = YpredictX_Q1_rho, ymax = YpredictX_Q3_rho), fill = "grey80", alpha = 0.5) +
    geom_ribbon(aes(ymin = XpredictY_Q1_rho, ymax = XpredictY_Q3_rho), fill = "grey80", alpha = 0.5) +
    geom_line(aes(y = YpredictX_rho, linetype = paste(y, "xmap", x))) +
    geom_line(aes(y = XpredictY_rho, linetype = paste(x, "xmap", y))) +
    geom_point(aes(y = YpredictX_rho)) +
    geom_point(aes(y = XpredictY_rho)) +

    labs(linetype = "Cross Map", subtitle = paste("Optimal Embedding Dimension of ", y, ":", Ey)) +
    theme_bw() +
    scale_x_continuous(breaks = seq(0, 600, by = 100)) +
    theme(text = element_text(size = 20), legend.position = c(.65, .65), plot.subtitle = element_text(size = 12)) +
    xlab("Library Size") +
    ylab("Forecaset Skill (ρ)")
  ggsave(paste0(ccm_rho_fig_path, x, " & ", y, "_BestE(", y, ")_", Ey, ".tiff"),
         units = "in",
         width = 8, height = 8, dpi = 300, compression = 'lzw')

  ggplot(sum_ccm_res_Ey, aes(libsize)) +
    geom_ribbon(aes(ymin = YpredictX_Q1_rmse, ymax = YpredictX_Q3_rmse), fill = "grey80", alpha = 0.5) +
    geom_ribbon(aes(ymin = XpredictY_Q1_rmse, ymax = XpredictY_Q3_rmse), fill = "grey80", alpha = 0.5) +
    geom_line(aes(y = YpredictX_rmse, linetype = paste(y, "xmap", x))) +
    geom_line(aes(y = XpredictY_rmse, linetype = paste(x, "xmap", y))) +
    geom_point(aes(y = YpredictX_rmse)) +
    geom_point(aes(y = XpredictY_rmse)) +

    labs(linetype = "Cross Map", subtitle = paste("Optimal Embedding Dimension of ", y, ":", Ey)) +
    theme_bw() +
    scale_x_continuous(breaks = seq(0, 600, by = 100)) +
    theme(text = element_text(size = 20), legend.position = c(.65, .65), plot.subtitle = element_text(size = 12)) +
    xlab("Library Size") +
    ylab("Forecaset Skill (RMSE)")
  ggsave(paste0(ccm_rmse_fig_path, x, " & ", y, "_BestE(", y, ")_", Ey, ".tiff"),
         units = "in",
         width = 8, height = 8, dpi = 300, compression = 'lzw')
}

ccm_colname <- c("A", "B")
ccm_func(ccm_colname)
# ccm_colname <- c("B", "A")
# ccm_func(ccm_colname)
ccm_colname <- c("AH1pdm", "AH3")
ccm_func(ccm_colname)
# ccm_colname <- c("AH3", "AH1pdm")
# ccm_func(ccm_colname)
ccm_colname <- c("Byama", "Bvic")
ccm_func(ccm_colname)
# ccm_colname <- c("Bvic", "Byama")
# ccm_func(ccm_colname)
ccm_colname <- c("AH1pdm", "B")
ccm_func(ccm_colname)
ccm_colname <- c("AH3", "B")
ccm_func(ccm_colname)
# ccm_block <- dt[, ccm_colname]

# for (i in 2:(length(dt)-1)) {
#   for (j in i + 1:(length(dt)-i)) {
#     # print(paste0(i,j))
#     x <- colnames(dt)[i]
#     y <- colnames(dt)[j]
#     print(paste(x,y))
#     ccm_block <- cbind(dt[x], dt[y])
#     ccm_func(ccm_block, x, y)
#
#   }
# }

ccm_colname <- c("AH1pdm", "AH3")
ccm_block <- dt[, ccm_colname]

x <- colnames(ccm_block[1])
# predicted col
y <- colnames(ccm_block[2])

# use Ex when use x to predict y
Ex <- sum_simp_res[which.min(sum_simp_res[[x]]),]$E
ccm_res_xtoy_Ex <- ccm(ccm_block, E = Ex, lib = c(1, NROW(ccm_block)), lib_column = 1, target_column = 2, lib_sizes = seq(Ex + 1, nrow(dt), by = 10))
ccm_res_tau_xtoy_Ex <- ccm(ccm_block, E = Ex, lib = c(1, NROW(ccm_block)), tau = -4:4,lib_column = 1, target_column = 2, lib_sizes = NROW(ccm_block))
ccm_stats_xtoy_Ex <- ccm_stats_func(ccm_res_xtoy_Ex)

# use Ey when use y to predict x
Ey <- sum_simp_res[which.min(sum_simp_res[[y]]),]$E
ccm_res_ytox_Ey <- ccm(ccm_block, E = Ey, lib = c(1, NROW(ccm_block)), lib_column = 2, target_column = 1, lib_sizes = seq(Ey + 1, nrow(dt), by = 10))
ccm_stats_ytox_Ey <- ccm_stats_func(ccm_res_ytox_Ey)

sum_ccm_res <- data.frame(
  libsize_Ex = ccm_stats_xtoy_Ex$lib_size,
  libsize_Ey = ccm_stats_ytox_Ey$lib_size,

  YpredictX_rho = ccm_stats_ytox_Ey$rho,
  YpredictX_rmse = ccm_stats_ytox_Ey$rmse,
  YpredictX_Q1_rho = ccm_stats_ytox_Ey$rhoQ1,
  YpredictX_Q3_rho = ccm_stats_ytox_Ey$rhoQ3,
  YpredictX_Q1_rmse = ccm_stats_ytox_Ey$rmseQ1,
  YpredictX_Q3_rmse = ccm_stats_ytox_Ey$rmseQ3,

  XpredictY_rho = ccm_stats_xtoy_Ex$rho,
  XpredictY_rmse = ccm_stats_xtoy_Ex$rmse,
  XpredictY_Q1_rho = ccm_stats_xtoy_Ex$rhoQ1,
  XpredictY_Q3_rho = ccm_stats_xtoy_Ex$rhoQ3,
  XpredictY_Q1_rmse = ccm_stats_xtoy_Ex$rmseQ1,
  XpredictY_Q3_rmse = ccm_stats_xtoy_Ex$rmseQ3
)

ggplot(sum_ccm_res) +
  geom_ribbon(aes(x= libsize_Ey,ymin = YpredictX_Q1_rho, ymax = YpredictX_Q3_rho), fill = "grey80", alpha = 0.5) +
  geom_ribbon(aes(x= libsize_Ex,ymin = XpredictY_Q1_rho, ymax = XpredictY_Q3_rho), fill = "grey80", alpha = 0.5) +
  geom_line(aes(x= libsize_Ey,y = YpredictX_rho, linetype = paste(y, "xmap", x))) +
  geom_line(aes(x= libsize_Ex,y = XpredictY_rho, linetype = paste(x, "xmap", y))) +
  geom_point(aes(x= libsize_Ey,y = YpredictX_rho)) +
  geom_point(aes(x= libsize_Ex,y = XpredictY_rho)) +

  labs(linetype = "Cross Map") +
  # theme_bw() +
  scale_x_continuous(breaks = seq(0, 600, by = 100)) +
  theme(text = element_text(size = 20), legend.position = c(.65, .65), plot.subtitle = element_text(size = 12)) +
  xlab("Library Size") +
  ylab("Cross Map Skill (ρ)")
ggsave(paste0(ccm_rho_fig_path, x, " & ", y, ".tiff"),
       units = "in",
       width = 8, height = 8, dpi = 300, compression = 'lzw')

ggplot(sum_ccm_res, aes(libsize)) +
  geom_ribbon(aes(ymin = YpredictX_Q1_rmse, ymax = YpredictX_Q3_rmse), fill = "grey80", alpha = 0.5) +
  geom_ribbon(aes(ymin = XpredictY_Q1_rmse, ymax = XpredictY_Q3_rmse), fill = "grey80", alpha = 0.5) +
  geom_line(aes(y = YpredictX_rmse, linetype = paste(y, "xmap", x))) +
  geom_line(aes(y = XpredictY_rmse, linetype = paste(x, "xmap", y))) +
  geom_point(aes(y = YpredictX_rmse)) +
  geom_point(aes(y = XpredictY_rmse)) +

  labs(linetype = "Cross Map", subtitle = paste("Optimal Embedding Dimension of ", x, ":", Ex)) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 600, by = 100)) +
  theme(text = element_text(size = 20), legend.position = c(.65, .65), plot.subtitle = element_text(size = 12)) +
  xlab("Library Size") +
  ylab("Forecaset Skill (RMSE)")
ggsave(paste0(ccm_rmse_fig_path, x, " & ", y, "_BestE(", x, ")_", Ex, ".tiff"),
       units = "in",
       width = 8, height = 8, dpi = 300, compression = 'lzw')

ggplot(sum_ccm_res_Ey, aes(libsize)) +
  geom_ribbon(aes(ymin = YpredictX_Q1_rho, ymax = YpredictX_Q3_rho), fill = "grey80", alpha = 0.5) +
  geom_ribbon(aes(ymin = XpredictY_Q1_rho, ymax = XpredictY_Q3_rho), fill = "grey80", alpha = 0.5) +
  geom_line(aes(y = YpredictX_rho, linetype = paste(y, "xmap", x))) +
  geom_line(aes(y = XpredictY_rho, linetype = paste(x, "xmap", y))) +
  geom_point(aes(y = YpredictX_rho)) +
  geom_point(aes(y = XpredictY_rho)) +

  labs(linetype = "Cross Map", subtitle = paste("Optimal Embedding Dimension of ", y, ":", Ey)) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 600, by = 100)) +
  theme(text = element_text(size = 20), legend.position = c(.65, .65), plot.subtitle = element_text(size = 12)) +
  xlab("Library Size") +
  ylab("Forecaset Skill (ρ)")
ggsave(paste0(ccm_rho_fig_path, x, " & ", y, "_BestE(", y, ")_", Ey, ".tiff"),
       units = "in",
       width = 8, height = 8, dpi = 300, compression = 'lzw')

ggplot(sum_ccm_res_Ey, aes(libsize)) +
  geom_ribbon(aes(ymin = YpredictX_Q1_rmse, ymax = YpredictX_Q3_rmse), fill = "grey80", alpha = 0.5) +
  geom_ribbon(aes(ymin = XpredictY_Q1_rmse, ymax = XpredictY_Q3_rmse), fill = "grey80", alpha = 0.5) +
  geom_line(aes(y = YpredictX_rmse, linetype = paste(y, "xmap", x))) +
  geom_line(aes(y = XpredictY_rmse, linetype = paste(x, "xmap", y))) +
  geom_point(aes(y = YpredictX_rmse)) +
  geom_point(aes(y = XpredictY_rmse)) +

  labs(linetype = "Cross Map", subtitle = paste("Optimal Embedding Dimension of ", y, ":", Ey)) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 600, by = 100)) +
  theme(text = element_text(size = 20), legend.position = c(.65, .65), plot.subtitle = element_text(size = 12)) +
  xlab("Library Size") +
  ylab("Forecaset Skill (RMSE)")
ggsave(paste0(ccm_rmse_fig_path, x, " & ", y, "_BestE(", y, ")_", Ey, ".tiff"),
       units = "in",
       width = 8, height = 8, dpi = 300, compression = 'lzw')

#
# ccm_colname <- c("B", "A")
#
# ccm_colname <- c("AH1pdm", "AH3")
#
# ccm_colname <- c("AH3", "AH1pdm")
#
# ccm_colname <- c("Byama", "Bvic")
#
# ccm_colname <- c("Bvic", "Byama")
#
# ccm_colname <- c("A", "B")
#
# ccm_block <- dt[, ccm_colname]
# BestE <- sum_simp_res[which.min(sum_simp_res[[colnames(ccm_block)[1]]]),]$E
# ccm_res_col1to2 <- ccm(ccm_block, E = BestE, lib = c(1, NROW(ccm_block)),
#                        lib_column = 1, target_column = 2,
#                        num_samples = 100,
#                        lib_sizes = seq(BestE + 1, nrow(dt), by = 20),
#                        stats_only = FALSE)
#
# ccm_means <- ccm_means(ccm_res_col1to2)
# quantile_res <- aggregate(rho ~ lib_size, data = ccm_res_col1to2, FUN = "quantile", probs = c(25, 75) / 100)
# ccm_quantile <- data.frame(quantile_res[[2]])
# ccm_stats <- cbind(ccm_means, ccm_quantile)
#
# sum_ccm_res <- data.frame(
#   libsize = ccm_stats$lib_size,
#   rho = ccm_stats$rho,
#   rmse = ccm_stats$rmse,
#   Q1 = ccm_stats$X25.,
#   Q3 = ccm_stats$X75.
# )
#
# ccm_plt <- ggplot(sum_ccm_res, aes(x = libsize)) +
#   geom_line(aes(y = rmse, colour = paste(colnames(ccm_block)[1],"xmap",colnames(ccm_block)[2])))

# ccm_colname <- c("Bvic", "Byama")
# ccm_block <- dt[, ccm_colname]
# write.csv(ccm_block,"/Users/yining/Downloads/ccm_testdata_B.csv")
#
# # x <- colnames(ccm_block[1])
# # # predicted col
# # y <- colnames(ccm_block[2])
# y1_ <- list()
# # use Ex when use x to predict y
# Ex <- sum_simp_res[which.min(sum_simp_res[[x]]),]$E
# x_sur <- rEDM::make_surrogate_data(ccm_block, method = "random_shuffle", num_surr = 100)
# y_sur <- rEDM::make_surrogate_data(ccm_block[2], method = "random_shuffle", num_surr = 100)
# for (i in 1:5){
#   ccm_block_surr <- cbind(ccm_block[x], y_sur[i])
#   ccm_res_xtoySurr <- ccm(ccm_block_surr, E = Ex, lib = c(1, NROW(ccm_block)), lib_column = 1, target_column = 2, lib_sizes = seq(Ex + 1, nrow(ccm_block), by = 20))
#   ccm_stats_xtoySurr_Ex <- rEDM::ccm_means(ccm_res_xtoySurr)
#   y1_[[i]] <- pmax(0, ccm_stats_xtoySurr_Ex$rmse)
# }


#
# ccm_res_xtoy_Ex_surrMin <- rEDM::ccm(ccm_block, E = Ex, lib = c(1, NROW(ccm_block)), target_column = 2,lib_sizes = Ex+1,num_samples=1000)
# ccm_res_xtoy_Ex_surrMax <- rEDM::ccm(ccm_block, E = Ex, lib = c(1, NROW(ccm_block)), target_column = 2,lib_sizes = nrow(ccm_block),num_samples=1000)
# quantile(ccm_res_xtoy_Ex_surrMin$rmse,p=c(0.025,0.5,0.975))
# quantile(ccm_res_xtoy_Ex_surrMax$rmse,p=c(0.025,0.5,0.975))


#
# simp_res <- simplex(y)
# Ey <- simp_res$E[which.min(simp_res$rmse)]
# plot(simp_res$E, simp_res$rmse, las = 1, xlab = "E", ylab = "RMSE")
# mtext(glue("The optimal embedding dim is {Ey}"), side = 3)
