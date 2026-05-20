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

#####################################################################
#################### INSTALL AND IMPORT PACKAGES ####################
#####################################################################
rm(list = ls())
# install.packages("remotes")
library(remotes)
#install.packages("pacman")
#install.packages("devtools")
#remotes::install_github("ha0ye/rEDM")
#remotes::install_github("yutakaos/rUIC")
#remotes::install_github("ong8181/macamts")

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
  ggdensity,
  paletteer
)
# Load library
packageVersion("rEDM") # v 0.7.5
# packageVersion("macam") # v 0.1.4
packageVersion("macamts") # v 0.1.4
packageVersion("rUIC") # v 0.9.12
# library(rEDM)
# vignette("rEDM-tutorial")
theme_set(theme_cowplot())

######################################################
##################### FUNCTIONS ######################
######################################################

# 有意なUIC結果をまとめて出す関数
process_df <- function(df, cause_var, effect_var) {
  res_95 <- subset(df, pval < 0.05 & ete > quantile_95)
  
  if (nrow(res_95) > 0) {
    signif_res <- res_95
    chosen_quantile <- "quantile_95"
  } else {
    res_90 <- subset(df, pval < 0.05 & ete > quantile_90)
    signif_res <- res_90
    chosen_quantile <- "quantile_90"
  }
  
  if (nrow(signif_res) == 0) {
    return(NULL)
  }
  
  selected_idx <- which.max(signif_res$ete)
  selected_row <- signif_res[selected_idx,]
  
  if (chosen_quantile == "quantile_95") {
    summary_signif_res_one <- data.frame(
      cause = cause_var,
      effected = effect_var,
      E = selected_row$E,
      tp = selected_row$tp,
      ete = selected_row$ete,
      quantile_90 = selected_row$quantile_90,
      quantile_95 = selected_row$quantile_95,
      weight = 3,  # pval > quantile_95
      stringsAsFactors = FALSE
    )
  } else {
    summary_signif_res_one <- data.frame(
      cause = cause_var,
      effected = effect_var,
      E = selected_row$E,
      tp = selected_row$tp,
      ete = selected_row$ete,
      quantile_90 = selected_row$quantile_90,
      quantile_95 = selected_row$quantile_95,
      weight = 1,  # pval > quantile_90
      stringsAsFactors = FALSE
    )
  }
  
  return(summary_signif_res_one)
}

# ransom_seedを入れたmake_surrogate_seasonal()関数
make_surrogate_seasonal_randomseed <- function(ts,
                                               num_surr = 100,
                                               T_period = 12,
                                               random_seed = NULL) {
  if (is.data.frame(ts)) {
    ts <- ts[[1]]
  }
  
  if (any(!is.finite(ts)))
    stop("input time series contained invalid values")
  
  n <- length(ts)
  I_season <- suppressWarnings(matrix(1:T_period, nrow = n, ncol = 1))
  
  # Calculate seasonal cycle using smooth.spline
  seasonal_F <- smooth.spline(c(I_season - T_period, I_season, I_season + T_period),
                              c(ts, ts, ts))
  seasonal_cyc <- predict(seasonal_F, I_season)$y
  seasonal_resid <- ts - seasonal_cyc
  
  # Set random seed if provided
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }
  
  # Generate surrogate data
  matrix(unlist(lapply(seq(num_surr), function(i) {
    seasonal_cyc + sample(seasonal_resid, n)
  })), ncol = num_surr)
}

# UICからS-mapまで解析を行う関数
# effect_var: 影響されたvar
# cause_var: 影響を与えるvar
# A cause Bの場合、A: cause_var, B: effect_var
UIC_Smap_func <- function(df, effect_var, cause_var, numSurr) {
  # df <- df_log
  # effect_var <- "flu"
  # cause_var <- "rs"
  # numSurr <- 3
  
  print(
    "---------------------------------------------------------------------------------------"
  )
  print(
    paste(
      "--------------- Effected:",
      effect_var,
      "| Cause:",
      cause_var,
      " ---------------"
    )
  )
  if (dir.exists(paste0(uic_surr_tbl_path,cause_var,"_cause_",effect_var,"_surrogate_data.csv")) == TRUE){
    return(NULL)
  }
  print(summary_signif_res) # 今まで有意なUIC結果を出力
  print("UIC analysis --- start")
  # print(head(df)) # dataframeを確認
  tp_range <- -12:0 # 時間のラグ（-tp）
  E_range <- 0:20 # 埋め込み次元（E）
  
  ## Does cause_var cause effect_var?
  uic_res <- uic.optimal(
    df,
    lib_var = effect_var,
    tar_var = cause_var,
    E = E_range,
    tau = 1,
    tp = tp_range,
    alpha = 0.05
  )
  BestE <- uic_res$E[1] + 1
  print(paste("Optimal E for", effect_var, ":", BestE))
  
  # save UIC results
  write.csv(
    uic_res,
    paste0(
      uic_surr_tbl_path,
      "result/",
      cause_var,
      "_cause_",
      effect_var,
      "_UIC_result.csv"
    )
  )
  print("save uic result --- done")
  print(
    paste0(
      "path : ",
      uic_surr_tbl_path,
      "result/",
      cause_var,
      "_cause_",
      effect_var,
      "_UIC_result.csv"
    )
  )
  
  ## Plot UIC results
  ggplot(uic_res, aes(x = tp, y = ete)) +
    geom_line() +
    labs(title = paste0("UIC (", cause_var, " causes ", effect_var, "?)")) +
    theme(title = element_text(face = "bold"), text = element_text(size = 20)) +
    geom_point(aes(color = pval < 0.05), size = 4) +
    theme(legend.position = "none") +
    scale_color_manual(values = c("black", "red")) +
    labs(x = "Time Lag (tp)", y = "Effective Transfer Entropy", color = "p < 0.05")
  # save plot
  ggsave(
    paste0(uic_res_path, cause_var, "_cause_", effect_var, "_UIC.tiff"),
    units = "in",
    width = 10,
    height = 8,
    dpi = 300,
    compression = 'lzw'
  )
  print("plot and save uic figure --- done")
  print(paste0(
    "path : ",
    uic_res_path,
    cause_var,
    "_cause_",
    effect_var,
    "_UIC.tiff"
  ))
  
  UIC_significant <- subset(uic_res, pval < 0.05)
  
  if (nrow(UIC_significant) == 0) {
    cat("significant UIC result: null\n")
    return(NULL)
  } else {
    cat("significant UIC result:\n")
    print(paste("Continue seasonality test using surrogate data"))
  }
  
  ## Generate seasonal surrogate data for effect_var
  # effect_surr <- rEDM::make_surrogate_seasonal(df[[effect_var]], num_surr = numSurr)
  effect_surr <- make_surrogate_seasonal_randomseed(df[[effect_var]], num_surr = numSurr, random_seed = 1234)
  # save surrogate data
  write.csv(
    effect_surr,
    paste0(
      uic_surr_tbl_path,
      cause_var,
      "_cause_",
      effect_var,
      "_surrogate_data.csv"
    )
  )
  print(
    paste0(
      "path : ",
      uic_surr_tbl_path,
      cause_var,
      "_cause_",
      effect_var,
      "_surrogate_data.csv"
    )
  )
  print("generate and save seasonal surrogate data for effect_var --- done")
  
  ## Compute UIC for surrogate data
  ete_surr <- data.frame(tp = tp_range)
  for (i in 1:ncol(effect_surr)) {
    block_tmp <- data.frame(effect = effect_surr[, i], cause = df[[cause_var]])
    ete_surr_i <- uic.optimal(
      block_tmp,
      lib_var = "effect",
      tar_var = "cause",
      E = E_range,
      tau = 1,
      tp = tp_range,
      num_surr = 1
    )
    res_ete_tp <- data.frame(ete_surr_i %>% select(tp, ete))
    colnames(res_ete_tp) <- c("tp", paste0("ete_", i))
    ete_surr <- merge(ete_surr, res_ete_tp, by = "tp")
  }
  # save
  write.csv(
    ete_surr,
    paste0(
      uic_surr_tbl_path,
      "result/",
      cause_var,
      "_cause_",
      effect_var,
      "_ete_tp_surrogate_data.csv"
    )
  )
  print("compute UIC result for surrogate data and save --- done")
  print(
    paste0(
      "path : ",
      uic_surr_tbl_path,
      "result/",
      cause_var,
      "_cause_",
      effect_var,
      "_ete_tp_surrogate_data.csv"
    )
  )
  
  ## Calculate quantiles from surrogate data
  quantile_list <- list()
  for (i in 1:nrow(ete_surr)) {
    ete_surr_noTP <- subset(ete_surr, select = -c(tp))
    res_quantile <- quantile(unlist(ete_surr_noTP[i, ]),
                             probs = c(0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99))
    quantile_list[[paste0("tp_", i - length(tp_range))]] <- res_quantile
  }
  df_quantiles <- data.frame(do.call(rbind, quantile_list))
  # save
  write.csv(
    df_quantiles,
    paste0(
      uic_surr_tbl_path,
      "result/",
      cause_var,
      "_cause_",
      effect_var,
      "_quantile_surrogate_data.csv"
    )
  )
  print("calculate quantiles from surrogate data and save --- done")
  print(
    paste0(
      "path : ",
      uic_surr_tbl_path,
      "result/",
      cause_var,
      "_cause_",
      effect_var,
      "_quantile_surrogate_data.csv"
    )
  )
  
  ## Determine significant UIC values (use the dataframe above)
  significant_UIC <- data.frame(
    E = uic_res$E,
    tp = uic_res$tp,
    ete = uic_res$ete,
    pval = uic_res$pval,
    quantile_90 = df_quantiles$X90.,
    quantile_95 = df_quantiles$X95.,
    quantile_97.5 = df_quantiles$X97.5.,
    quantile_99 = df_quantiles$X99.
  )
  write.csv(
    significant_UIC,
    paste0(
      uic_surr_tbl_path,
      cause_var,
      "_cause_",
      effect_var,
      "_uic_p_surr_res.csv"
    )
  )
  # print(paste("quantiles df: tp =", significant_UIC$tp,
  #             "; ete =", significant_UIC$ete,
  #             "; quantile 95 =", significant_UIC$quantile_95))
  print(
    paste0(
      "saving --- done; path : ",
      uic_surr_tbl_path,
      cause_var,
      "_cause_",
      effect_var,
      "_uic_p_surr_res.csv"
    )
  )
  
  ## Plot significant UIC results
  ggplot(significant_UIC, aes(x = tp)) +
    geom_line(aes(y = ete), linetype = "solid") +
    geom_line(aes(y = quantile_95), linetype = "longdash") +
    geom_line(aes(y = quantile_90), linetype = "dotted") +
    labs(title = paste0("UIC (", cause_var, " causes ", effect_var, "?)")) +
    theme(title = element_text(face = "bold"), text = element_text(size = 20)) +
    geom_point(aes(y = ete, color = pval < 0.05), size = 4) +
    scale_x_continuous(breaks = seq(-12, 0, by = 4)) +
    theme(legend.position = "none") +
    scale_color_manual(values = c("black", "red")) +
    labs(x = "Time Lag (tp)", y = "Effective Transfer Entropy", color = "p < 0.05")
  
  ggsave(
    paste0(
      uicsurr_res_path,
      cause_var,
      "_cause_",
      effect_var,
      "_UIC_surr.tiff"
    ),
    units = "in",
    width = 10,
    height = 8,
    dpi = 300,
    compression = 'lzw'
  )
  print("plot and save uic_surrogate figure --- done")
  print(
    paste0(
      "path : ",
      uicsurr_res_path,
      cause_var,
      "_cause_",
      effect_var,
      "_UIC_surr.tiff"
    )
  )
  
  # S-map
  uic_surr_significant_res <- read.csv(
    paste0(
      uic_surr_tbl_path,
      cause_var,
      "_cause_",
      effect_var,
      "_uic_p_surr_res.csv"
    )
  )
  uic_res_pval <- uic_surr_significant_res[uic_surr_significant_res$pval < 0.05, ]
  
  tmp_res <- process_df(uic_res_pval, cause_var, effect_var)
  
  if (is.null(tmp_res)) {
    cat("significant result: null\n")
    return(NULL)
  } else {
    cat("significant result:\n")
    BestTP <- tmp_res$tp
    print(paste("Best tp:", BestTP))
    print(tmp_res)
  }
  
  print("S-map analysis --- start")
  # Create lagged variables
  if (BestTP == 0) {
    lagged_effect <- dplyr::lead(df[[effect_var]], n = 1) # tp = 0
  } else {
    lagged_effect <- dplyr::lag(df[[effect_var]], n = abs(BestTP) - 1)
  }
  lagged_cause <- make_block(df[[cause_var]], max_lag = BestE)
  
  # Combine into a state-space block
  smap_block <- cbind(lagged_effect, lagged_cause[, 2:ncol(lagged_cause)])
  write.csv(
    smap_block,
    paste0(
      res_save_path,
      "table/smap/block/",
      effect_var,
      "_effected_by_",
      cause_var,
      ".csv"
    )
  )
  
  ## Optimize the theta parameter
  print("S-map parameter search --- start")
  theta_range <- c(0,
                   1e-04,
                   3e-04,
                   0.001,
                   0.003,
                   0.01,
                   0.03,
                   0.1,
                   0.3,
                   0.5,
                   0.75,
                   1,
                   1.5,
                   2,
                   3,
                   4,
                   6,
                   8)
  lambda_range <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
  # theta_range <- c(0, 3e-04,0.003, 0.1, 0.5, 1, 3, 8)
  # lambda_range <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1)
  stat_res <- data.frame(
    N = numeric(),
    theta = numeric(),
    lambda = numeric(),
    rho = numeric(),
    mae = numeric(),
    rmse = numeric(),
    stringsAsFactors = FALSE
  )
  
  # for (theta in theta_range) {
  #   rsmap_ridge <- macamts::extended_lnlp(smap_block, theta = theta, lambda = 0.1,
  #                                         regularized = TRUE, alpha = 0, random_seed = 1234, save_smap_coefficients = TRUE)
  #   stat_res <- rbind(stat_res, rsmap_ridge$stats)
  # }
  
  for (theta_value in theta_range) {
    # print(paste0("theta: ",theta_value))
    for (lambda_value in lambda_range) {
      # print(paste0("lambda: ",lambda_value))
      rsmap_ridge <- macamts::extended_lnlp(
        smap_block,
        theta = theta_value,
        lambda = lambda_value,
        regularized = TRUE,
        alpha = 0,
        random_seed = 1234,
        save_smap_coefficients = TRUE
      )
      stat_res_one <- data.frame(
        N = rsmap_ridge$stats$N,
        theta = theta_value,
        lambda = lambda_value,
        rho = rsmap_ridge$stats$rho,
        mae = rsmap_ridge$stats$mae,
        rmse = rsmap_ridge$stats$rmse,
        stringsAsFactors = FALSE
      )
      stat_res <- rbind(stat_res, stat_res_one)
    }
  }
  
  BestTheta <- stat_res[which.min(stat_res$rmse), "theta"]
  BestLambda <- stat_res[which.min(stat_res$rmse), "lambda"]
  tmp_res$Theta <- BestTheta
  tmp_res$Lambda <- BestLambda
  print(
    paste(
      "Optimal theta&labmda for",
      effect_var,
      "and",
      cause_var,
      ":",
      BestTheta,
      ";",
      BestLambda
    )
  )
  write.csv(
    stat_res,
    paste0(
      res_save_path,
      "table/smap/parameter/",
      effect_var,
      "_effected_by_",
      cause_var,
      "_thetaLambda",
      ".csv"
    )
  )
  
  ## Perform regularized S-map analysis with optimized parameters
  rsmap_ridge <- macamts::extended_lnlp(
    smap_block,
    theta = BestTheta,
    lambda = BestLambda,
    regularized = TRUE,
    alpha = 0,
    random_seed = 1234,
    save_smap_coefficients = TRUE
  )
  smap_pred_res <- rsmap_ridge$model_output
  write.csv(
    smap_pred_res,
    paste0(
      res_save_path,
      "table/smap/pred_res/",
      effect_var,
      "_effected_by_",
      cause_var,
      ".csv"
    )
  )
  
  coef_res <- rsmap_ridge$smap_coefficients
  write.csv(
    coef_res,
    paste0(
      res_save_path,
      "table/smap/coef/",
      effect_var,
      "_effected_by_",
      cause_var,
      ".csv"
    )
  )
  
  ## Plot observed vs. predicted values
  maxValue <- max(max(na.omit(smap_pred_res$obs)), max(na.omit(smap_pred_res$pred))) + 1
  ggplot(smap_pred_res, aes(x = obs, y = pred, color = time)) +
    geom_abline(slope = 1,
                linetype = "dashed",
                color = "black") +
    geom_point() +
    xlim(c(NA, maxValue)) +
    ylim(c(NA, maxValue)) +
    labs(x = "Observed", y = "Predicted", color = "Time")
  
  ggsave(
    paste0(
      xmapping_fig_path,
      "pred_obs_",
      cause_var,
      "_causes_",
      effect_var,
      ".tiff"
    ),
    units = "in",
    width = 8,
    height = 8,
    dpi = 300,
    compression = 'lzw'
  )
  
  ## Prepare coefficients for plotting
  coef_res_cols <- data.frame(time = coef_res$time,
                              c_effect = coef_res$c_1,
                              c_cause = coef_res[[paste0("c_", BestE + 1)]])
  write.csv(
    coef_res_cols,
    paste0(
      res_save_path,
      "table/smap/coef/",
      "summary_",
      effect_var,
      "_effected_by_",
      cause_var,
      ".csv"
    )
  )
  
  coef_average <- mean(coef_res_cols$c_cause, na.rm = TRUE)
  tmp_res$coef_average <- coef_average
  coef_median <- median(coef_res_cols$c_cause, na.rm = TRUE)
  tmp_res$coef_median <- coef_median
  coef_variance <- var(coef_res_cols$c_cause, na.rm = TRUE)
  tmp_res$coef_variance <- coef_variance
  coef_std_dev <- sd(coef_res_cols$c_cause, na.rm = TRUE)
  tmp_res$coef_std_dev <- coef_std_dev
  
  if (coef_average > 0) {
    tmp_res$color <- "red"
  } else if (coef_average == 0) {
    tmp_res$color <- "grey"
  } else {
    tmp_res$color <- "blue"
  }
  
  
  ## Create mathematical expressions for plotting
  c_effect_expr <- substitute(frac(partialdiff * x1, partialdiff * x2),
                              list(x1 = bquote(.(
                                as.symbol(effect_var)
                              )(t + 1)), x2 = bquote(.(
                                as.symbol(effect_var)
                              )(t))))
  c_cause_expr <- substitute(frac(partialdiff * x1, partialdiff * x2),
                             list(x1 = bquote(.(
                               as.symbol(effect_var)
                             )(t + 1)), x2 = bquote(.(
                               as.symbol(cause_var)
                             )(t))))
  
  ## Reshape data for plotting
  coef_res_long <- gather(coef_res_cols,
                          key = "coef",
                          value = "value",
                          c_effect,
                          c_cause)
  
  ## Plot coefficients
  ggplot(coef_res_long, aes(x = coef, y = value)) +
    geom_boxplot() +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.5) +
    scale_x_discrete(labels = c("c_effect" = c_effect_expr, "c_cause" = c_cause_expr)) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "gray") +
    theme(
      text = element_text(size = 20),
      legend.position = "none",
      axis.title.x = element_blank()
    ) +
    labs(y = "Regularized S-map Coefficient")
  
  ggsave(
    paste0(
      fig_path,
      "smapFig/",
      "coefficients_",
      cause_var,
      "_causes_",
      effect_var,
      ".tiff"
    ),
    units = "in",
    width = 8,
    height = 8,
    dpi = 300,
    compression = 'lzw'
  )
  
  ## Time series plots
  date_xaxis <- as.POSIXct(df$Date, format = "%Y-%m-%d")
  
  # Observed and predicted plot
  ObsPred_plt <- ggplot(smap_pred_res, aes(date_xaxis)) +
    geom_line(aes(y = obs, linetype = "Observation")) +
    geom_line(aes(y = pred, linetype = "Prediction")) +
    scale_colour_grey(start = 0, end = 0.8) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    ylab("Incidence") +
    labs(linetype = effect_var) +
    theme(
      text = element_text(size = 20),
      axis.title.x = element_blank(),
      legend.position = c(0.80, 0.9),
      axis.title.y = element_text(angle = 0, vjust = 0.5)
    )
  
  ggsave(
    paste0(
      fig_path,
      "smapFig/",
      "ObsPred_",
      cause_var,
      "_causes_",
      effect_var,
      ".tiff"
    ),
    units = "in",
    width = 24,
    height = 8,
    dpi = 300,
    compression = 'lzw'
  )
  
  # Coefficient plots over time
  c_effect_plot <- ggplot(coef_res_cols, aes(date_xaxis)) +
    geom_line(aes(y = c_effect)) +
    geom_hline(
      yintercept = quantile(na.omit(coef_res_cols$c_effect), c(0.25, 0.5, 0.75)),
      linetype = "dashed",
      color = c("gray", "black", "gray")
    ) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    labs(y = c_effect_expr) +
    theme(
      text = element_text(size = 20),
      axis.title.x = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5)
    )
  
  ggsave(
    paste0(
      fig_path,
      "smapFig/",
      "effect_",
      cause_var,
      "_causes_",
      effect_var,
      ".tiff"
    ),
    units = "in",
    width = 24,
    height = 8,
    dpi = 300,
    compression = 'lzw'
  )
  
  c_cause_plot <- ggplot(coef_res_cols, aes(date_xaxis)) +
    geom_line(aes(y = c_cause)) +
    geom_hline(
      yintercept = quantile(na.omit(coef_res_cols$c_cause), c(0.25, 0.5, 0.75)),
      linetype = "dashed",
      color = c("gray", "black", "gray")
    ) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    labs(y = c_cause_expr) +
    theme(
      text = element_text(size = 20),
      axis.title.x = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5)
    )
  
  ggsave(
    paste0(
      fig_path,
      "smapFig/",
      "cause_",
      cause_var,
      "_causes_",
      effect_var,
      ".tiff"
    ),
    units = "in",
    width = 24,
    height = 8,
    dpi = 300,
    compression = 'lzw'
  )
  
  # Combine plots if needed
  merged_plot <- plot_grid(ObsPred_plt,
                           c_effect_plot,
                           c_cause_plot,
                           align = "v",
                           ncol = 1)
  # merged_plot
  ggsave(
    paste0(
      fig_path,
      "smapFig/",
      "merged_",
      cause_var,
      "_causes_",
      effect_var,
      ".tiff"
    ),
    units = "in",
    width = 24,
    height = 24,
    dpi = 300,
    compression = 'lzw'
  )
  
  print(tmp_res)
  return(tmp_res)
}

#####################################################################
##################### SET PATH AND IMPORT DATA ######################
#####################################################################
getwd()

## Set path
dir.create(
  "result/oknw_ID/oknw_ID_weather_result/"
)
global_res_save_path <- paste0(
  "result/oknw_ID/oknw_ID_weather_result/",
  Sys.Date()
)
dir.create(file.path(global_res_save_path))
print(global_res_save_path) # check the path

dt_ID <- read.csv(
  "data/oknw/oknwID_10to19_mainIsland_2islands_250207.csv"
)
dt_ID <- dt_ID %>%
  mutate(date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1")))
dt_ID$date <- as.Date(dt_ID$date) # date is for time series visualization
colnames(dt_ID)

regions <- c("Miyako","Yaeyama", "mainIsland")
for (region_name in regions) {
  region <- region_name
  #region <- "mainIsland"
  
  # log(): base "e"
  #
  
  dt0 <- data.frame(
    Date = dt_ID$date,
    flu = dt_ID[[paste0("flu_", region)]],
    temp = dt_ID[[paste0("meanTemp_C_", region)]],
    humid = dt_ID[[paste0("meanHumidity_percent_", region)]],
    rainfall = dt_ID[[paste0("SumRain_mm_", region)]]
  )
  
  df_log <- data.frame(
    Date = dt_ID$date,
    flu = log(dt_ID[[paste0("flu_", region)]] + 1),
    temp = log(dt_ID[[paste0("meanTemp_C_", region)]] + 1),
    humid = log(dt_ID[[paste0("meanHumidity_percent_", region)]] + 1),
    rainfall = log(dt_ID[[paste0("SumRain_mm_", region)]] + 1)
  )
  
  summary_res <- data.frame(
    colname = character(),
    length = numeric(),
    zero_count = numeric(),
    summary_log = numeric(),
    summary_raw = numeric(),
    mean = numeric(),
    
    min = numeric(),
    p25 = numeric(),
    median = numeric(),
    p75 = numeric(),
    max = numeric()
  )
  for (cols in colnames(df_log)[2:length(df_log)]) {
    print(cols)
    
    summary_res_one <- data.frame(
      colname = cols,
      length = length(df_log[[cols]]),
      zero_count = sum(df_log[[cols]] == 0),
      summary_log = sum(df_log[[cols]]),
      mean = mean(df_log[[cols]]),
      min = min(df_log[[cols]]),
      p25 = quantile(df_log[[cols]], 0.25),
      median = median(df_log[[cols]]),
      p75 = quantile(df_log[[cols]], 0.75),
      max = max(df_log[[cols]])
    )
    print(summary_res_one)
    summary_res <- rbind(summary_res, summary_res_one)
  }
  print(summary_res)
  write.csv(summary_res,paste0(global_res_save_path, "dataframe_summary_log_",region,".csv"))
  
  summary_res <- data.frame(
    colname = character(),
    length = numeric(),
    zero_count = numeric(),
    summary_log = numeric(),
    summary_raw = numeric(),
    mean = numeric(),
    
    min = numeric(),
    p25 = numeric(),
    median = numeric(),
    p75 = numeric(),
    max = numeric()
  )
  for (cols in colnames(dt0)[2:length(dt0)]) {
    print(cols)
    
    summary_res_one <- data.frame(
      colname = cols,
      length = length(dt0[[cols]]),
      zero_count = sum(dt0[[cols]] == 0),
      summary_raw = sum(dt0[[cols]]),
      mean = mean(dt0[[cols]]),
      min = min(dt0[[cols]]),
      p25 = quantile(dt0[[cols]], 0.25),
      median = median(dt0[[cols]]),
      p75 = quantile(dt0[[cols]], 0.75),
      max = max(dt0[[cols]])
    )
    print(summary_res_one)
    summary_res <- rbind(summary_res, summary_res_one)
  }
  print(summary_res)
  #
  write.csv(summary_res,paste0(global_res_save_path, "dataframe_summary_raw_",region,".csv"))
  
  # データをlong形式に変換
  data_long <- pivot_longer(dt0, cols = c("flu", "temp", "humid", "rainfall"),
                            names_to = "variable", values_to = "value")
  
  # パネルの順序を設定（上からflu, temp, humid, rainfall）
  data_long$variable <- factor(data_long$variable, levels = c("flu", "temp", "humid", "rainfall"))
  
  # ggplotでDateをx軸、valueをy軸として、変数ごとにパネル分けします
  ggplot(data_long, aes(x = Date, y = value)) +
    geom_line() +
    facet_wrap(~ variable, ncol = 1, scales = "free_y") +
    scale_x_date(date_breaks = "3 month", date_labels = "%Y-%m") +
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1), 
          axis.text.y = element_text(size = 24)) +
    labs(x = "Date")
  ggsave(paste0("result/oknw_ID/oknw_ID_weather_result/raw_ts/",region,"_raw_ts_panel.tiff"),
         units = "in", width = 38, height = 28, dpi = 300, compression = 'lzw')

  ggplot(dt0, aes(Date)) +
    geom_point(aes(y = flu, colour = "influenza")) +
    geom_point(aes(y = temp, colour = "temprature")) +
    geom_point(aes(y = humid, colour = "humidity")) +
    geom_point(aes(y = rainfall, colour = "precipitation ")) +
    
    geom_line(aes(y = flu, colour = "influenza")) +
    geom_line(aes(y = temp, colour = "temprature")) +
    geom_line(aes(y = humid, colour = "humidity")) +
    geom_line(aes(y = rainfall, colour = "precipitation ")) +
    
    labs(colour = "Influenza") +
    scale_x_date(date_breaks = "6 month", date_labels = "%Y-%m") +
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size = 24, angle = 45, hjust = 1), 
          axis.text.y = element_text(size = 24)) +
    scale_color_paletteer_d("fishualize::Acanthurus_sohal") +
    xlab("Date") +
    ylab("Incidence")
  ggsave(paste0("result/oknw_ID/oknw_ID_weather_result/raw_ts/",region,"_ts.tiff"),
         units = "in", width = 38, height = 7, dpi = 300, compression = 'lzw')
  
  ggplot(df_log, aes(Date)) +
    geom_point(aes(y = flu, colour = "influenza")) +
    geom_point(aes(y = temp, colour = "temprature")) +
    geom_point(aes(y = humid, colour = "humidity")) +
    geom_point(aes(y = rainfall, colour = "precipitation ")) +
    
    geom_line(aes(y = flu, colour = "influenza")) +
    geom_line(aes(y = temp, colour = "temprature")) +
    geom_line(aes(y = humid, colour = "humidity")) +
    geom_line(aes(y = rainfall, colour = "precipitation ")) +
    
    labs(colour = "Influenza") +
    scale_x_date(date_breaks = "6 month", date_labels = "%Y-%m") +
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size = 24, angle = 45, hjust = 1), 
          axis.text.y = element_text(size = 24)) +
    scale_color_paletteer_d("fishualize::Acanthurus_sohal") +
    xlab("Date") +
    ylab("Incidence")
  ggsave(paste0("result/oknw_ID/oknw_ID_weather_result/raw_ts/",region,"_log_ts.tiff"),
         units = "in", width = 38, height = 7, dpi = 300, compression = 'lzw')
  
  ggplot(dt0, aes(Date)) +
    geom_point(aes(y = flu, colour = "influenza")) +
    geom_line(aes(y = flu, colour = "influenza")) +
    labs(colour = "Influenza") +
    scale_x_date(date_breaks = "6 month", date_labels = "%Y-%m") +
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size = 24, angle = 45, hjust = 1), 
          axis.text.y = element_text(size = 24)) +
    xlab("Date") +
    ylab("Incidence")
  ggsave(paste0("result/oknw_ID/oknw_ID_weather_result/raw_ts/",region,"_flu_raw_ts.tiff"),
         units = "in", width = 38, height = 7, dpi = 300, compression = 'lzw')
  
  ggplot(dt0, aes(Date)) +
    geom_point(aes(y = temp, colour = "temprature")) +
    geom_line(aes(y = temp, colour = "temprature")) +
    
    scale_x_date(date_breaks = "6 month", date_labels = "%Y-%m") +
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size = 24, angle = 45, hjust = 1), 
          axis.text.y = element_text(size = 24)) +
    xlab("Date") +
    ylab("Incidence")
  ggsave(paste0("result/oknw_ID/oknw_ID_weather_result/raw_ts/",region,"_temp_raw_ts.tiff"),
         units = "in", width = 38, height = 7, dpi = 300, compression = 'lzw')
  
  ggplot(dt0, aes(Date)) +
    geom_point(aes(y = humid, colour = "humidity")) +
    geom_line(aes(y = humid, colour = "humidity")) +
    
    scale_x_date(date_breaks = "6 month", date_labels = "%Y-%m") +
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size = 24, angle = 45, hjust = 1), 
          axis.text.y = element_text(size = 24)) +
    xlab("Date") +
    ylab("Incidence")
  ggsave(paste0("result/oknw_ID/oknw_ID_weather_result/raw_ts/",region,"_humid_raw_ts.tiff"),
         units = "in", width = 38, height = 7, dpi = 300, compression = 'lzw')
  
  ggplot(dt0, aes(Date)) +
    geom_point(aes(y = rainfall, colour = "precipitation")) +
    geom_line(aes(y = rainfall, colour = "precipitation")) +
    
    scale_x_date(date_breaks = "6 month", date_labels = "%Y-%m") +
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size = 24, angle = 45, hjust = 1), 
          axis.text.y = element_text(size = 24)) +
    xlab("Date") +
    ylab("Incidence")
  ggsave(paste0("result/oknw_ID/oknw_ID_weather_result/raw_ts/",region,"_rain_raw_ts.tiff"),
         units = "in", width = 38, height = 7, dpi = 300, compression = 'lzw')
  
  
}

# 並列計算のプランを設定 (必要に応じて変更)
# plan(multisession, workers = future::availableCores() - 1)
#print(paste("available core:", availableCores() - 1))

# regions <- c("Miyako", "mainIsland", "Yaeyama")
# regions <- c("Miyako","Yaeyama", "mainIsland")
regions <- c("mainIsland")

for (region_name in regions) {
  region <- region_name
  # region <- "Miyako"
  
  ID_name <- c(
    "flu",
    "rs",
    "pcf",
    "gas",
    "aig",
    "cpox",
    "hfmd",
    "eryt",
    "es",
    "herp",
    "mumps",
    "ahc",
    "ekc",
    "bm",
    "am",
    "mp",
    "cp"
  )
  df_log <- data.frame(
    Date = dt_ID$date,
    flu = log(dt_ID[[paste0("flu_", region)]] + 1),
    rs = log(dt_ID[[paste0("rs_", region)]] + 1),
    pcf = log(dt_ID[[paste0("pcf_", region)]] + 1),
    gas = log(dt_ID[[paste0("gas_", region)]] + 1),
    aig = log(dt_ID[[paste0("aig_", region)]] + 1),
    cpox = log(dt_ID[[paste0("cpox_", region)]] + 1),
    hfmd = log(dt_ID[[paste0("hfmd_", region)]] + 1),
    eryt = log(dt_ID[[paste0("eryt_", region)]] + 1),
    es = log(dt_ID[[paste0("es_", region)]] + 1),
    herp = log(dt_ID[[paste0("herp_", region)]] + 1),
    mumps = log(dt_ID[[paste0("mumps_", region)]] + 1),
    ahc = log(dt_ID[[paste0("ahc_", region)]] + 1),
    ekc = log(dt_ID[[paste0("ekc_", region)]] + 1),
    bm = log(dt_ID[[paste0("bm_", region)]] + 1),
    am = log(dt_ID[[paste0("am_", region)]] + 1),
    mp = log(dt_ID[[paste0("mp_", region)]] + 1),
    cp = log(dt_ID[[paste0("cp_", region)]] + 1),
    temp = log(dt_ID[[paste0("meanTemp_C_", region)]] + 1),
    humid = log(dt_ID[[paste0("meanHumidity_percent_", region)]] + 1),
    rainfall = log(dt_ID[[paste0("SumRain_mm_", region)]] + 1)
  )
  
  for (cols in colnames(df_log)[2:length(df_log)]) {
    len <- length(df_log[[cols]])
    sum_var <- sum(df_log[[cols]])
    sum_zero <- sum(df_log[[cols]] == 0)
    print(paste(
      "colname:",
      cols,
      "; length:",
      len,
      "; zero count:",
      sum_zero
    ))
    if (sum_var == 0) {
      df_log <- df_log[, colnames(df_log) != cols]
    }
    if (sum_zero > len * 0.8) {
      df_log <- df_log[, colnames(df_log) != cols]
    }
    if (ncol(df_log) < 5) {
      return(NULL)
    }
  }
  #
  print(colnames(df_log))
  
  area <- region
  print(paste0("------------------------------ ",area," ------------------------------"))
  ## global save path
  res_save_path <- paste0(global_res_save_path, "/", area, "/")
  dir.create(file.path(res_save_path))
  print(res_save_path)
  print(dir.exists(res_save_path))
  ## figure save path
  fig_path <- paste0(res_save_path, "figure/")
  dir.create(file.path(fig_path))
  print(fig_path)
  print(dir.exists(fig_path))
  # UIC figure result save path
  xmapping_fig_path <- paste0(fig_path, "crossMapping/")
  dir.create(file.path(xmapping_fig_path))
  print(xmapping_fig_path)
  print(dir.exists(xmapping_fig_path))
  #
  uic_fig_path <- paste0(fig_path, "uicFig/")
  dir.create(file.path(uic_fig_path))
  print(uic_fig_path)
  print(dir.exists(uic_fig_path))
  #
  uic_res_path <- paste0(uic_fig_path, "uic/")
  dir.create(file.path(uic_res_path))
  print(uic_res_path)
  print(dir.exists(uic_res_path))
  #
  uicsurr_res_path <- paste0(uic_fig_path, "uic_surr/")
  dir.create(file.path(uicsurr_res_path))
  print(uicsurr_res_path)
  print(dir.exists(uicsurr_res_path))
  # sMap figure result save path
  smap_fig_path <- paste0(fig_path, "smapFig/")
  dir.create(file.path(smap_fig_path))
  print(smap_fig_path)
  print(dir.exists(smap_fig_path))
  ## table save path
  tbl_path <- paste0(res_save_path, "table/")
  dir.create(file.path(tbl_path))
  print(tbl_path)
  print(dir.exists(tbl_path))
  # uic table save path
  uic_tbl_path <- paste0(tbl_path, "uic/")
  dir.create(file.path(uic_tbl_path))
  print(uic_tbl_path)
  print(dir.exists(uic_tbl_path))
  #
  uic_res_tbl_path <- paste0(uic_tbl_path, "result/")
  dir.create(file.path(uic_res_tbl_path))
  print(uic_res_tbl_path)
  print(dir.exists(uic_res_tbl_path))
  #
  uic_surr_tbl_path <- paste0(uic_tbl_path, "surrogate_dt/")
  dir.create(file.path(uic_surr_tbl_path))
  print(uic_surr_tbl_path)
  print(dir.exists(uic_surr_tbl_path))
  dir.create(file.path(paste0(uic_surr_tbl_path, "result/")))
  print(dir.exists(paste0(uic_surr_tbl_path, "result/")))
  # sMap table save path
  smap_tbl_path <- paste0(tbl_path, "smap/")
  dir.create(file.path(smap_tbl_path))
  print(smap_tbl_path)
  print(dir.exists(smap_tbl_path))
  #
  smap_parameter_tbl_path <- paste0(smap_tbl_path, "parameter/")
  dir.create(file.path(smap_parameter_tbl_path))
  print(smap_parameter_tbl_path)
  print(dir.exists(smap_parameter_tbl_path))
  #
  smap_block_tbl_path <- paste0(smap_tbl_path, "block/")
  dir.create(file.path(smap_block_tbl_path))
  print(smap_block_tbl_path)
  print(dir.exists(smap_block_tbl_path))
  #
  smap_coef_tbl_path <- paste0(smap_tbl_path, "coef/")
  dir.create(file.path(smap_coef_tbl_path))
  print(smap_coef_tbl_path)
  print(dir.exists(smap_coef_tbl_path))
  #
  smap_predres_tbl_path <- paste0(smap_tbl_path, "pred_res/")
  dir.create(file.path(smap_predres_tbl_path))
  print(smap_predres_tbl_path)
  print(dir.exists(smap_predres_tbl_path))
  
  # 解析結果を保存するdataframe：結果はpythonのnetworkxパッケージに直接入力したら、相互作用のネットワークが作成できる
  # dataframeのイメージ：（実際、"tp, weight, cause, effected, color"がnetworkxに入力として使われる）
  summary_signif_res <- data.frame(
    cause = character(),
    # cause_var
    effected = character(),
    # effect_var: 影響されるvar
    E = numeric(),
    # Best E
    tp = numeric(),
    # Best tp: 有意な結果の中に、eteが一番高いtpの値
    ete = numeric(),
    # Best tpの時のeteの値
    quantile_90 = numeric(),
    # サローゲートデータの結果の90% quantile (p = 0.1)
    quantile_95 = numeric(),
    # サローゲートデータの結果の95% quantile (p = 0.05)
    weight = numeric(),
    # p < 0.1の時、weightは1；p < 0.05の時、weightは3
    Theta = numeric(),
    # BestTheta from S-map
    Lambda = numeric(),
    # BestLambda from S-map
    color = character(),
    # s-mapのcoefficientの結果の平均がpositive: red; negative: blue
    coef_average = numeric(),
    # s-mapのcoefficientの結果の平均
    coef_median = numeric(),
    # s-mapのcoefficientの結果の中央値
    coef_variance = numeric(),
    # s-mapのcoefficientの結果の分散
    coef_std_dev = numeric(),
    # s-mapのcoefficientの結果の標準偏差
    stringsAsFactors = FALSE
  )
  
  set.seed(1234)
    
  weather <- c("temp", "humid", "rainfall")
  col_names <-  colnames(df_log)
  ID_names <- col_names[!col_names %in% c("Date", "temp", "humid", "rainfall")]
  comb <- expand.grid(
    effect_var = ID_names,
    cause_var = weather,
    stringsAsFactors = FALSE
  )
  comb <- comb[comb$effect_var != comb$cause_var, ]
  for (i in row(comb)) {
    effect_var <- comb$effect_var[i]
    cause_var  <- comb$cause_var[i]
    numSurr    <- 2000
    temp_res <- UIC_Smap_func(df_log, effect_var, cause_var, numSurr)
    if (!is.null(temp_res)) {
      summary_signif_res <- rbind(summary_signif_res, temp_res)
    }
    print(summary_signif_res)
    write.csv(summary_signif_res, paste0(res_save_path, "table/", "significant_res_summary_temp.csv"))
    
    }
      
    write.csv(
      summary_signif_res,
      file = paste0(
        res_save_path,
        "table/",
        "significant_res_summary_final.csv"
      ),
      row.names = FALSE
    )
}   
