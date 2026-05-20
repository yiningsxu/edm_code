rm(list = ls())
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

library(remotes)
remotes::install_github("ha0ye/rEDM")
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
  writexl
)
packageVersion("rEDM")
# library(rEDM)
# vignette("rEDM-tutorial")


## Set path
# setwd handled by edm_code bootstrap
## plot save path
fig_date_path <- "result/figure/oknw_ID/24-09-11/"
# # Raw time series
# raw_ts_path <- paste0(fig_date_path , "incidence_oknw_ID.tiff")
# # Standalized time series
# standalized_ts_path <- paste0(fig_date_path ,"standarlized_incidence_oknw_ID.tiff")
# # Best E
# bestE_fig_path <- paste0(fig_date_path ,"BestE_H1H3B.tiff")
# # nonlinear theta
# theta_fig_path <- paste0(fig_date_path ,"theta_H1H3B.tiff")
# attractor
attractor_path <- paste0(fig_date_path, "ID_attractor/")
# interaction boxplot
boxplot_path <- paste0(fig_date_path, "ID_interaction/boxplot/")
# interaction time series
ts_path <- paste0(fig_date_path, "ID_interaction/ts/")
# ccm rho
ccm_rho_fig_path <- paste0(fig_date_path, "ID_causality/rho_")
# ccm rmse
ccm_rmse_fig_path <- paste0(fig_date_path, "ID_causality/rmse_")
## --------------------- Preparation --------------------- ##
# ## import data
# dt0 <- read.csv("data/oknw/merged_incidence_weather_plt.csv")
# # str(dt0)
# # dt <- subset(dt, select = -c(X))
# dt0$date <- as.Date(dt0$date)
# # dt0 <- dt0[, c(147, 3:146)]
# # dt0
#
#
# dt <- data.frame(dt0[, 7:length(dt0)])
# dt_clean <- dt[, apply(dt, 2, function(x) !any(is.na(x)))]
# dt_stad <- data.frame(scale(dt_clean))
# dt_stad <- mutate(dt_stad, Date = dt0$date)
# dt_stad <- dt_stad[, c(length(dt_stad), 1:length(dt_stad) - 1)]
# write.csv(dt_stad, "data/oknw/standarlized_oknw_ID_10to19.csv")

data <- read.csv("data/oknw/standarlized_oknw_ID_10to19.csv")
# 519 x 155
data <- subset(data, select = -c(X))
data$Date <- as.Date(data$Date)
str(data)
colnames(data)
data <- data %>% select_if(negate(anyNA))
# 468 x 133

# naha <- mutate(dt[, grep("_naha$", names(dt), value = TRUE)], Date = dt$Date)
# naha_dt <- naha[, c(length(naha), 1:length(naha) - 1)]
#
# north <- mutate(dt[, grep("_north$", names(dt), value = TRUE)], Date = data$Date)
# north_dt <- north[, c(length(north), 1:length(north) - 1)]
#
# middle <- mutate(dt[, grep("_middle$", names(dt), value = TRUE)], Date = dt$Date)
# middle_dt <- middle[, c(length(middle), 1:length(middle) - 1)]
#
# south <- mutate(dt[, grep("_south$", names(dt), value = TRUE)], Date = dt$Date)
# south_dt <- south[, c(length(south), 1:length(south) - 1)]
#
# myk <- mutate(dt[, grep("_myk$", names(dt), value = TRUE)], Date = dt$Date)
# myk_dt <- myk[, c(length(myk), 1:length(myk) - 1)]
#
# yeym <- mutate(dt[, grep("_yeym$", names(dt), value = TRUE)], Date = dt$Date)
# yeym_dt <- yeym[, c(length(yeym), 1:length(yeym) - 1)]
#
# oknw <- mutate(dt[, grep("_oknw$", names(dt), value = TRUE)], Date = dt$Date)
# oknw_dt <- oknw[, c(length(oknw), 1:length(oknw) - 1)]
#
# jp <- mutate(dt[, grep("_jp$", names(dt), value = TRUE)], Date = dt$Date)
# jp_dt <- jp[, c(length(jp), 1:length(jp) - 1)]

## --------------------- Multi S-map Result Table --------------------- ##
# multi_smap_res <- function(smap_block, x, y, i) {
#   multiSmap_res <- block_lnlp(smap_block, method = "s-map",
#                               theta = seq(0, 15, 0.5),
#                               target_column = x,
#                               first_column_time = FALSE,
#                               silent = TRUE, stats_only = FALSE)
#   multiSmap_min_theta <- multiSmap_res[which.min(multiSmap_res$rmse), "theta"]
#   multiSmap_result <- block_lnlp(smap_block, method = "s-map",
#                                  theta = multiSmap_min_theta,
#                                  target_column = x,
#                                  first_column_time = FALSE,
#                                  silent = TRUE, save_smap_coefficients = TRUE)
#
#   coef_res <- data.frame(multiSmap_result$smap_coefficients[[1]])
#   print(paste0(x, "(Y) & ", y, " (X)"))
#   multiSmap_res_dt[i, "counter"] <<- i
#   multiSmap_res_dt[i, "X"] <<- paste0(y)
#   multiSmap_res_dt[i, "Y"] <<- paste0(x)
#   multiSmap_res_dt[i, "Y_mean"] <<- mean(coef_res$c_1, na.rm = TRUE)
#   multiSmap_res_dt[i, "X_mean"] <<- mean(coef_res$c_2, na.rm = TRUE)
#   multiSmap_res_dt[i, "Y_median"] <<- median(coef_res$c_1, na.rm = TRUE)
#   multiSmap_res_dt[i, "X_median"] <<- median(coef_res$c_2, na.rm = TRUE)
#   multiSmap_res_dt[i, "Y_max"] <<- max(coef_res$c_1, na.rm = TRUE)
#   multiSmap_res_dt[i, "X_max"] <<- max(coef_res$c_2, na.rm = TRUE)
#   multiSmap_res_dt[i, "Y_min"] <<- min(coef_res$c_1, na.rm = TRUE)
#   multiSmap_res_dt[i, "X_min"] <<- min(coef_res$c_2, na.rm = TRUE)
# }

col_name <- c("naha_dt", "north_dt", "middle_dt", "south_dt", "myk_dt", "yeym_dt", "oknw_dt", "jp_dt")
col_pattern <- c("_Naha$", "_North$", "_Middle$", "_South$", "_Miyako$", "_Yaeyama$", "_Okinawa$", "_Japan$")
for (num in 1:length(col_pattern)) {
  dt0 <- mutate(data[, grep(col_pattern[num], names(data), value = TRUE)], Date = data$Date)
  print(col_pattern[num])
  print(paste(nrow(dt0),"x",ncol(dt0)))
  print(colnames(dt0))
}

for (num in 1:length(col_pattern)) {
  num <- num
  dt0 <- mutate(data[, grep(col_pattern[num], names(data), value = TRUE)], Date = data$Date)
  dt <- dt0[, c(length(dt0), 1:length(dt0) - 1)]
  mulSmap_resPath <- paste0("result/analysis_table/", col_name[num], "_multiSmap_res_dt.xlsx")

  simplex_res_bestE <- data.frame(E = simplex(dt[[colnames(dt)[2]]], E = 1:20, silent = T)$E)
  for (i in 2:length(dt)) {
    simplex_res <- simplex(dt[[colnames(dt)[i]]],
                           E = 1:20,
                           silent = T,
                           stats_only = FALSE)$rmse
    simplex_res_bestE[[colnames(dt)[i]]] <- simplex_res
  }
  sum_simp_res <- simplex_res_bestE

  simplex_res_list <- list()
  for (i in 2:length(dt)) {
    E_value <- sum_simp_res[which.min(sum_simp_res[[colnames(dt)[i]]]), "E"]
    simplex_res <- simplex(dt[[colnames(dt)[i]]],
                           E = E_value,
                           stats_only = FALSE)
    simplex_res_list[[colnames(dt)[i]]] <- list(simplex_res)
  }

  multiSmap_Xmean <- data.frame()
  multiSmap_Xmedian <- data.frame()
  multiSmap_Ymean <- data.frame()
  multiSmap_Ymedian <- data.frame()

  multiSmap_XmeanPath <- paste0("result/analysis_table/subpref/Xmean/", col_name[num], "_multiSmap_Xmean.xlsx")
  multiSmap_XmedianPath <- paste0("result/analysis_table/subpref/Xmedian/", col_name[num], "_multiSmap_Xmedian.xlsx")
  multiSmap_YmeanPath <- paste0("result/analysis_table/subpref/Ymean/", col_name[num], "_multiSmap_Ymean.xlsx")
  multiSmap_YmedianPath <- paste0("result/analysis_table/subpref/Ymedian/", col_name[num], "_multiSmap_Ymedian.xlsx")

  counter <- 1
  for (i in 2:length(dt)) {
    for (j in 2:length(dt)) {
      print(paste("counter:", counter))
      x <- colnames(dt)[j]
      y <- colnames(dt)[i]
      smap_block <- cbind(dt[x], dt[y])
      multiSmap_res <- block_lnlp(smap_block, method = "s-map",
                                  theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5,0.75, 1, 1.5, 2, 3, 4,5, 6, 7,8),
                                  target_column = 2,
                                  first_column_time = FALSE,
                                  silent = TRUE, stats_only = FALSE)
      multiSmap_min_theta <- multiSmap_res[which.min(multiSmap_res$rmse), "theta"]
      multiSmap_result <- block_lnlp(smap_block, method = "s-map",
                                     theta = multiSmap_min_theta,
                                     target_column = 2,
                                     first_column_time = FALSE,
                                     silent = TRUE, save_smap_coefficients = TRUE)

      coef_res <- data.frame(multiSmap_result$smap_coefficients[[1]])

      multiSmap_Xmean[j - 1, "X"] <- x
      multiSmap_Xmedian[j - 1, "X"] <- x
      multiSmap_Ymean[j - 1, "X"] <- x
      multiSmap_Ymedian[j - 1, "X"] <- x

      multiSmap_Xmean[j - 1, y] <- mean(coef_res$c_1, na.rm = TRUE)
      multiSmap_Xmedian[j - 1, y] <- median(coef_res$c_1, na.rm = TRUE)
      multiSmap_Ymean[j - 1, y] <- mean(coef_res$c_2, na.rm = TRUE)
      multiSmap_Ymedian[j - 1, y] <- median(coef_res$c_2, na.rm = TRUE)

      counter <- counter + 1
    }
  }
  write_xlsx(multiSmap_Xmean, path = multiSmap_XmeanPath)
  write_xlsx(multiSmap_Xmedian, path = multiSmap_XmedianPath)
  write_xlsx(multiSmap_Ymean, path = multiSmap_YmeanPath)
  write_xlsx(multiSmap_Ymedian, path = multiSmap_YmedianPath)
}


## -------------------------------------------------------------------- ##
## -------------------------------------------------------------------- ##
## -------------------------------------------------------------------- ##

## ----------------------------- data ----------------------------- ##
col_name <- c("naha_dt", "north_dt", "middle_dt", "south_dt", "myk_dt", "yeym_dt", "oknw_dt", "jp_dt")
col_pattern <- c("_naha$", "_north$", "_middle$", "_south$", "_myk$", "_yeym$", "_oknw$", "_jp$")
num <- 1
dt0 <- mutate(data[, grep(col_pattern[num], names(data), value = TRUE)], Date = data$Date)
dt <- dt0[, c(length(dt0), 1:length(dt0) - 1)]
mulSmap_resPath <- paste0("result/analysis_table/", col_name[num], "_multiSmap_res_dt.xlsx")

## --------------------- Simple Projection --------------------- ##
simplex_res_bestE <- data.frame(E = simplex(dt[[colnames(dt)[2]]], E = 1:20, silent = T)$E)
for (i in 2:length(dt)) {
  simplex_res <- simplex(dt[[colnames(dt)[i]]],
                         E = 1:20,
                         silent = T,
                         stats_only = FALSE)$rmse
  simplex_res_bestE[[colnames(dt)[i]]] <- simplex_res
}
sum_simp_res <- simplex_res_bestE


## --------------------- Result of Simplex Projection --------------------- ##
simplex_res_list <- list()
for (i in 2:length(dt)) {
  E_value <- sum_simp_res[which.min(sum_simp_res[[colnames(dt)[i]]]), "E"]
  simplex_res <- simplex(dt[[colnames(dt)[i]]],
                         E = E_value,
                         stats_only = FALSE)
  simplex_res_list[[colnames(dt)[i]]] <- list(simplex_res)
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

sum_smap_res <- data.frame("theta" = data.frame(uniSmap_res_list[[colnames(dt)[2]]])$theta)
for (i in 2:ncol(dt)) {
  print(colnames(dt)[i])
  rmse_value <- data.frame(uniSmap_res_list[[colnames(dt)[i]]])$rmse
  rmse <- data.frame(rmse_value)
  colnames(rmse) <- colnames(dt)[i]
  sum_smap_res <- cbind(sum_smap_res, rmse)
}

sum_smap_res

### --------------------------- Multivariate S-map --------------------------- ##
# multi_smap_func <- function(smap_block, x, y) {
#   multiSmap_res <- block_lnlp(smap_block, method = "s-map",
#                               theta = seq(0, 15, 0.5),
#                               target_column = x,
#                               first_column_time = FALSE,
#                               silent = TRUE, stats_only = FALSE)
#   multiSmap_min_theta <- multiSmap_res[which.min(multiSmap_res$rmse), "theta"]
#   multiSmap_result <- block_lnlp(smap_block, method = "s-map",
#                                  theta = multiSmap_min_theta,
#                                  target_column = x,
#                                  first_column_time = FALSE,
#                                  silent = TRUE, save_smap_coefficients = TRUE)
#
#   coef_res <- data.frame(multiSmap_result$smap_coefficients[[1]])
#
#   coef_res_long <- gather(coef_res, key = "variable", value = "value", c_1, c_2)
#   boxplot <- ggplot(coef_res_long, aes(x = variable, y = value)) +
#     geom_boxplot() +
#     theme_bw() +
#     scale_x_discrete(labels = c(
#       "c_1" = expression(frac(partialdiff * Y[t + 1], partialdiff * Y[t])),
#       "c_2" = expression(frac(partialdiff * Y[t + 1], partialdiff * X[t]))
#     )) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
#     theme(text = element_text(size = 20), legend.position = c(.25, .65), axis.title.x = element_blank(),) +
#     labs(x = paste0("Interactions Between ", x, " and ", y), y = "S-map coefficients")
#   ggsave(paste0(boxplot_path, x, "(Y) & ", y, " (X)", ".tiff"),
#          plot = boxplot, units = "in",
#          width = 5, height = 8, dpi = 300, compression = 'lzw')
#
#   date_xaxis <- as.POSIXct(dt$Date, format = "%Y-%m-%d")
#   c1Y_plt <- ggplot(coef_res, aes(date_xaxis)) +
#     geom_line(aes(y = c_1, linetype = paste0(x))) +
#     # geom_line(aes(y = c_2, linetype = "c_2")) +
#     scale_colour_grey(start = 0, end = .8) +
#     scale_x_datetime(
#       date_breaks = "1 year",
#       date_labels = "%Y",
#       minor_breaks = "1 month"
#     ) +
#     theme_bw() +
#     labs(y = expression(frac(partialdiff * Y[t + 1], partialdiff * Y[t]))) +
#     theme(text = element_text(size = 20), axis.title.x = element_blank(),
#           axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
#
#   c2X_plt <- ggplot(coef_res, aes(date_xaxis)) +
#     geom_line(aes(y = c_2, linetype = paste0(y))) +
#     scale_colour_grey(start = 0, end = .8) +
#     labs(y = expression(frac(partialdiff * Y[t + 1], partialdiff * X[t]))) +
#     scale_x_datetime(
#       date_breaks = "1 year",
#       date_labels = "%Y",
#       minor_breaks = "1 month"
#     ) +
#     theme_bw() +
#     theme(text = element_text(size = 20), axis.title.x = element_blank(),
#           axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
#
#   ggsave(paste0(ts_path, x, "&", y, "_c1_", x, ".tiff"),
#          plot = c1Y_plt, units = "in",
#          width = 14, height = 4, dpi = 300, compression = 'lzw')
#   ggsave(paste0(ts_path, x, "&", y, "_c2_", y, ".tiff"),
#          plot = c2X_plt, units = "in",
#          width = 14, height = 4, dpi = 300, compression = 'lzw')
# }
#
# for (i in 2:length(dt)) {
#   for (j in 2:length(dt)) {
#     if (i != j) {
#       x <- colnames(dt)[i]
#       y <- colnames(dt)[j]
#       smap_block <- cbind(dt[x], dt[y])
#       multi_smap_func(smap_block, x, y)
#     }
#   }
# }

# --------------------------------------------------------------- #
multiSmap_Xmean <- data.frame()
multiSmap_Xmedian <- data.frame()
multiSmap_Ymean <- data.frame()
multiSmap_Ymedian <- data.frame()

multiSmap_XmeanPath <- paste0("result/analysis_table/", col_name[num], "_multiSmap_Xmean.xlsx")
multiSmap_XmedianPath <- paste0("result/analysis_table/", col_name[num], "_multiSmap_Xmedian.xlsx")
multiSmap_YmeanPath <- paste0("result/analysis_table/", col_name[num], "_multiSmap_Ymean.xlsx")
multiSmap_YmedianPath <- paste0("result/analysis_table/", col_name[num], "_multiSmap_Ymedian.xlsx")

counter <- 1
for (i in 2:length(dt)) {
  for (j in 2:length(dt)) {
    print(paste("counter:", counter))
    x <- colnames(dt)[j]
    y <- colnames(dt)[i]
    smap_block <- cbind(dt[x], dt[y])
    multiSmap_res <- block_lnlp(smap_block, method = "s-map",
                                theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5,0.75, 1, 1.5, 2, 3, 4,5, 6, 7,8),
                                target_column = 2,
                                first_column_time = FALSE,
                                silent = TRUE, stats_only = FALSE)
    multiSmap_min_theta <- multiSmap_res[which.min(multiSmap_res$rmse), "theta"]
    multiSmap_result <- block_lnlp(smap_block, method = "s-map",
                                   theta = multiSmap_min_theta,
                                   target_column = 2,
                                   first_column_time = FALSE,
                                   silent = TRUE, save_smap_coefficients = TRUE)

    coef_res <- data.frame(multiSmap_result$smap_coefficients[[1]])

    multiSmap_Xmean[j - 1, y] <- mean(coef_res$c_1, na.rm = TRUE)
    multiSmap_Xmedian[j - 1, y] <- median(coef_res$c_1, na.rm = TRUE)
    multiSmap_Ymean[j - 1, y] <- mean(coef_res$c_2, na.rm = TRUE)
    multiSmap_Ymedian[j - 1, y] <- median(coef_res$c_2, na.rm = TRUE)



    rownames(multiSmap_res_Xmean)[j - 1] <- x
    counter <- counter + 1
  }
}
write_xlsx(multiSmap_res_dt, path = mulSmap_resPath)


#
# matrix(vector(), ncol = length(dt)-1)
# colnames(multiSmap_res_dt) <- c("counter", "vars", "X_mean", "Y_mean", "X_median", "Y_median", "X_max", "X_min", "Y_max", "Y_min")

multiSmap_res <- block_lnlp(smap_block, method = "s-map",
                            # theta = seq(0, 15, 0.5),
                            target_column = y,
                            first_column_time = FALSE,
                            silent = TRUE, stats_only = FALSE)
multiSmap_min_theta <- multiSmap_res[which.min(multiSmap_res$rmse), "theta"]
multiSmap_result <- block_lnlp(smap_block, method = "s-map",
                               theta = multiSmap_min_theta,
                               target_column = y,
                               first_column_time = FALSE,
                               silent = TRUE, save_smap_coefficients = TRUE)

coef_res <- data.frame(multiSmap_result$smap_coefficients[[1]])
rownames(multiSmap_res_dt)[j - 1] <- x
multiSmap_res_dt[j - 1, y] <- mean(coef_res$c_1, na.rm = TRUE)

multi_smap_res <- function(smap_block, x, y) {
  multiSmap_res <- block_lnlp(smap_block, method = "s-map",
                              # theta = seq(0, 15, 0.5),
                              target_column = y,
                              first_column_time = FALSE,
                              silent = TRUE, stats_only = FALSE)
  multiSmap_min_theta <- multiSmap_res[which.min(multiSmap_res$rmse), "theta"]
  multiSmap_result <- block_lnlp(smap_block, method = "s-map",
                                 theta = multiSmap_min_theta,
                                 target_column = y,
                                 first_column_time = FALSE,
                                 silent = TRUE, save_smap_coefficients = TRUE)

  coef_res <- data.frame(multiSmap_result$smap_coefficients[[1]])
  # print(paste0(x, "(Y) & ", y, " (X)"))
  # multiSmap_res_dt[i, "counter"] <<- i
  # multiSmap_res_dt[i, "vars"] <<- paste0(x, "(Y) & ", y, " (X)")
  # multiSmap_res_dt[i, "Y_mean"] <<- mean(coef_res$c_1, na.rm = TRUE)
  # multiSmap_res_dt[i, "X_mean"] <<- mean(coef_res$c_2, na.rm = TRUE)
  # multiSmap_res_dt[i, "Y_median"] <<- median(coef_res$c_1, na.rm = TRUE)
  # multiSmap_res_dt[i, "X_median"] <<- median(coef_res$c_2, na.rm = TRUE)
  # multiSmap_res_dt[i, "Y_max"] <<- max(coef_res$c_1, na.rm = TRUE)
  # multiSmap_res_dt[i, "X_max"] <<- max(coef_res$c_2, na.rm = TRUE)
  # multiSmap_res_dt[i, "Y_min"] <<- min(coef_res$c_1, na.rm = TRUE)
  # multiSmap_res_dt[i, "X_min"] <<- min(coef_res$c_2, na.rm = TRUE)
}

counter <- 1
for (i in 2:length(dt)) {
  start_time <- proc.time()
  for (j in 2:length(dt)) {
    if (i != j) {
      print(paste("counter:", counter))
      x <- colnames(dt)[i]
      y <- colnames(dt)[j]
      smap_block <- cbind(dt[x], dt[y])
      multi_smap_res(smap_block, x, y, counter)
      counter <- counter + 1
    }
  }
  end_time <- proc.time()
  time_taken <- end_time - start_time
  print(time_taken)
}

write_xlsx(multiSmap_res_dt, path = mulSmap_resPath)

## --------------------------- CCM --------------------------- ##
ccm_func <- function(ccm_block, x, y) {

  Ex <- sum_simp_res[which.min(sum_simp_res[[x]]),]$E
  Ey <- sum_simp_res[which.min(sum_simp_res[[y]]),]$E

  ccm_res_xtoy_Ex <- ccm(ccm_block, E = Ex, lib = c(1, NROW(ccm_block)), target_column = 2, lib_sizes = seq(Ex + 1, nrow(dt), by = 20))
  ccm_means_xtoy_Ex <- ccm_means(ccm_res_xtoy_Ex)
  ccm_res_xtoy_Ey <- ccm(ccm_block, E = Ey, lib = c(1, NROW(ccm_block)), target_column = 2, lib_sizes = seq(Ey + 1, nrow(dt), by = 20))
  ccm_means_xtoy_Ey <- ccm_means(ccm_res_xtoy_Ey)

  ccm_res_ytox_Ex <- ccm(ccm_block, E = Ex, lib = c(1, NROW(ccm_block)), target_column = 1, lib_sizes = seq(Ex + 1, nrow(dt), by = 20))
  ccm_means_ytox_Ex <- ccm_means(ccm_res_ytox_Ex)
  ccm_res_ytox_Ey <- ccm(ccm_block, E = Ey, lib = c(1, NROW(ccm_block)), target_column = 1, lib_sizes = seq(Ey + 1, nrow(dt), by = 20))
  ccm_means_ytox_Ey <- ccm_means(ccm_res_ytox_Ey)

  sum_ccm_res_Ex <- data.frame(
    libsize = ccm_means_ytox_Ex$lib_size,
    YpredictX_rho = ccm_means_ytox_Ex$rho,
    YpredictX_rmse = ccm_means_ytox_Ex$rmse,
    XpredictY_rho = ccm_means_xtoy_Ex$rho,
    XpredictY_rmse = ccm_means_xtoy_Ex$rmse
  )

  sum_ccm_res_Ey <- data.frame(
    libsize = ccm_means_ytox_Ey$lib_size,
    YpredictX_rho = ccm_means_ytox_Ey$rho,
    YpredictX_rmse = ccm_means_ytox_Ey$rmse,
    XpredictY_rho = ccm_means_xtoy_Ey$rho,
    XpredictY_rmse = ccm_means_xtoy_Ey$rmse
  )

  ggplot(sum_ccm_res_Ex, aes(libsize)) +
    geom_line(aes(y = YpredictX_rho, linetype = paste(y, "predicts", x))) +
    geom_line(aes(y = XpredictY_rho, linetype = paste(x, "predicts", y))) +
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
    geom_line(aes(y = YpredictX_rmse, linetype = paste(y, "predicts", x))) +
    geom_line(aes(y = XpredictY_rmse, linetype = paste(x, "predicts", y))) +
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
    geom_line(aes(y = YpredictX_rho, linetype = paste(y, "predicts", x))) +
    geom_line(aes(y = XpredictY_rho, linetype = paste(x, "predicts", y))) +
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
    geom_line(aes(y = YpredictX_rmse, linetype = paste(y, "predicts", x))) +
    geom_line(aes(y = XpredictY_rmse, linetype = paste(x, "predicts", y))) +
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

x <- colnames(dt)[3]
y <- colnames(dt)[4]
ccm_block <- cbind(dt[x], dt[y])
ccm_func(ccm_block, x, y)

for (i in 2:(length(dt) - 1)) {
  for (j in i + 1:(length(dt) - i)) {
    # print(paste0(i,j))
    x <- colnames(dt)[i]
    y <- colnames(dt)[j]
    print(paste(x, y))
    ccm_block <- cbind(dt[x], dt[y])
    ccm_func(ccm_block, x, y)

  }
}