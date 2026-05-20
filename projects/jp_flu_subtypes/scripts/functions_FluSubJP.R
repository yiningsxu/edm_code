# Function to perform UIC analysis with surrogate data
Surr_UIC <- function(effect_var, cause_var, numSurr) {
  tp_range <- -12:0
  E_range <- 0:20

  ## Does cause_var cause effect_var?
  uic_res <- uic.optimal(df, lib_var = effect_var, tar_var = cause_var, E = E_range, tau = 1, tp = tp_range, alpha = 0.05)
  print(paste("Optimal E for", effect_var, ":", uic_res$E[1] + 1))

  # Save UIC results
  write.csv(uic_res, paste0(uic_surr_tbl_path, "result/", cause_var, "_cause_", effect_var, "_UIC_result.csv"))

  # Plot UIC results
  ggplot(uic_res, aes(x = tp, y = ete)) +
    geom_line() +
    labs(title = paste0("UIC (", cause_var, " causes ", effect_var, "?)")) +
    theme(title = element_text(face = "bold"), text = element_text(size = 20)) +
    geom_point(aes(color = pval < 0.05), size = 4) +
    theme(legend.position = "none") +
    scale_color_manual(values = c("black", "red")) +
    labs(x = "Time Lag (tp)", y = "Effective Transfer Entropy", color = "p < 0.05")

  ggsave(paste0(uic_plot_path, "uic/", cause_var, "_cause_", effect_var, "_UIC.tiff"),
         units = "in", width = 10, height = 8, dpi = 300, compression = 'lzw')

  ## Generate seasonal surrogate data for effect_var
  effect_surr <- rEDM::make_surrogate_seasonal(df[[effect_var]], num_surr = numSurr)
  write.csv(effect_surr, paste0(uic_surr_tbl_path, "surrogate_dt/", cause_var, "_cause_", effect_var, "_surrogate_data.csv"))

  ## Compute UIC for surrogate data
  ete_surr <- data.frame(tp = tp_range)
  for (i in 1:ncol(effect_surr)) {
    block_tmp <- data.frame(effect = effect_surr[, i], cause = df[[cause_var]])
    ete_surr_i <- uic.optimal(block_tmp,
                              lib_var = "effect",
                              tar_var = "cause", E = E_range, tau = 1, tp = tp_range, num_surr = 1)
    res_ete_tp <- data.frame(ete_surr_i %>% select(tp, ete))
    colnames(res_ete_tp) <- c("tp", paste0("ete_", i))
    ete_surr <- merge(ete_surr, res_ete_tp, by = "tp")
  }
  write.csv(ete_surr, paste0(uic_surr_tbl_path, "result/", cause_var, "_cause_", effect_var, "_ete_tp_surrogate_data.csv"))

  ## Calculate quantiles from surrogate data
  quantile_list <- list()
  for (i in 1:nrow(ete_surr)) {
    ete_surr_noTP <- subset(ete_surr, select = -c(tp))
    res_quantile <- quantile(unlist(ete_surr_noTP[i, ]), probs = c(0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99))
    quantile_list[[paste0("tp_", i - length(tp_range) + 1)]] <- res_quantile
  }
  df_quantiles <- data.frame(do.call(rbind, quantile_list))
  write.csv(df_quantiles, paste0(uic_surr_tbl_path, "result/", cause_var, "_cause_", effect_var, "_quantile_surrogate_data.csv"))

  ## Determine significant UIC values
  significant_UIC <- data.frame(
    tp = uic_res$tp,
    ete = uic_res$ete,
    pval = uic_res$pval,
    quantile_90 = df_quantiles$X90.,
    quantile_95 = df_quantiles$X95.,
    quantile_97.5 = df_quantiles$X97.5.,
    quantile_99 = df_quantiles$X99.
  )
  write.csv(significant_UIC, paste0(uic_surr_tbl_path, cause_var, "_cause_", effect_var, "_uic_p_surr_res.csv"))

  ## Plot significant UIC results
  ggplot(significant_UIC, aes(x = tp)) +
    geom_line(aes(y = ete), linetype = "solid") +
    geom_line(aes(y = quantile_95), linetype = "longdash") +
    geom_line(aes(y = quantile_90), linetype = "dotted") +
    labs(title = paste0("UIC (", cause_var, " causes ", effect_var, "?)")) +
    theme(title = element_text(face = "bold"), text = element_text(size = 20)) +
    geom_point(aes(y = ete, color = pval < 0.05), size = 4) +
    theme(legend.position = "none") +
    scale_color_manual(values = c("black", "red")) +
    labs(x = "Time Lag (tp)", y = "Effective Transfer Entropy", color = "p < 0.05")

  ggsave(paste0(uic_plot_path, "uic_surr/", cause_var, "_cause_", effect_var, "_UIC_surr.tiff"),
         units = "in", width = 10, height = 8, dpi = 300, compression = 'lzw')
}




# effect_var <- "A_H1N1"
# cause_var<- "B"
# Function to perform regularized S-map analysis
reg_smap_func <- function(effect_var, cause_var) {
  print(paste("Effected variable:", effect_var, "| Cause variable:", cause_var))

  ## Calculate the optimal embedding dimension (E)
  tp_range <- -12:0
  E_range <- 0:20
  uic_res <- uic.optimal(df, lib_var = effect_var, tar_var = cause_var, E = E_range, tau = 1, tp = tp_range, alpha = 0.05)
  BestE <- uic_res$E[1] + 1
  print(paste("Optimal E for", effect_var, ":", BestE))

  ## Prepare lagged variables based on the best time lag (tp)
  # uic_surr_tbl_path <- "result/FluSub_JP/2024-10-08/table/uic/"
  uic_res <- read.csv(paste0(uic_surr_tbl_path, cause_var, "_cause_", effect_var, "_uic_p_surr_res.csv"))
  uic_res_pval <- uic_res[uic_res$pval < 0.05, ]
  BestTP <- uic_res_pval[which.max(uic_res_pval$ete), ]$tp
  print(paste("Best tp:", BestTP))

  # Create lagged variables
  if (BestTP == 0) {
    lagged_effect <- dplyr::lead(df[[effect_var]], n = 1) # tp = 0
  } else {
    lagged_effect <- dplyr::lag(df[[effect_var]], n = abs(BestTP) - 1)
  }
  lagged_cause <- make_block(df[[cause_var]], max_lag = BestE)

  # Combine into a state-space block
  smap_block <- cbind(lagged_effect, lagged_cause[, 2:ncol(lagged_cause)])
  write.csv(smap_block, paste0(res_save_path, "table/smap/block/", effect_var, "_effected_by_", cause_var, ".csv"))

  ## Optimize the theta parameter
  theta_range <- c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
  stat_res <- list()
  for (theta in theta_range) {
    rsmap_ridge <- macamts::extended_lnlp(smap_block, theta = theta, lambda = 0.1,
                                          regularized = TRUE, alpha = 0, random_seed = 1234, save_smap_coefficients = TRUE)
    stat_res <- rbind(stat_res, rsmap_ridge$stats)
  }
  stat_res <- cbind(theta_range, stat_res)
  BestTheta <- stat_res[which.min(stat_res$rmse), "theta_range"]
  print(paste("Optimal theta for", effect_var, "and", cause_var, ":", BestTheta))

  ## Perform regularized S-map analysis with optimized parameters
  rsmap_ridge <- macamts::extended_lnlp(smap_block, theta = BestTheta, lambda = 0.1,
                                        regularized = TRUE, alpha = 0, random_seed = 1234, save_smap_coefficients = TRUE)
  smap_pred_res <- rsmap_ridge$model_output
  write.csv(smap_pred_res, paste0(res_save_path, "table/smap/pred_res/", effect_var, "_effected_by_", cause_var, ".csv"))

  coef_res <- rsmap_ridge$smap_coefficients
  write.csv(coef_res, paste0(res_save_path, "table/smap/coef/", effect_var, "_effected_by_", cause_var, ".csv"))

  ## Plot observed vs. predicted values
  maxValue <- max(max(na.omit(smap_pred_res$obs)), max(na.omit(smap_pred_res$pred))) + 1
  ggplot(smap_pred_res, aes(x = obs, y = pred, color = time)) +
    geom_abline(slope = 1, linetype = "dashed", color = "black") +
    geom_point() +
    xlim(c(NA, maxValue)) +
    ylim(c(NA, maxValue)) +
    labs(x = "Observed", y = "Predicted", color = "Time")

  ggsave(paste0(xmap_plot_path, "pred_obs_", cause_var, "_causes_", effect_var, ".tiff"),
         units = "in", width = 8, height = 8, dpi = 300, compression = 'lzw')

  ## Prepare coefficients for plotting
  coef_res_cols <- data.frame(
    time = coef_res$time,
    c_effect = coef_res$c_1,
    c_cause = coef_res[[paste0("c_", BestE + 1)]]
  )
  write.csv(coef_res_cols, paste0(res_save_path, "table/smap/coef/", "summary_", effect_var, "_effected_by_", cause_var, ".csv"))

  ## Create mathematical expressions for plotting
  c_effect_expr <- substitute(frac(partialdiff * x1, partialdiff * x2), list(x1 = bquote(.(as.symbol(effect_var))(t + 1)), x2 = bquote(.(as.symbol(effect_var))(t))))
  c_cause_expr <- substitute(frac(partialdiff * x1, partialdiff * x2), list(x1 = bquote(.(as.symbol(effect_var))(t + 1)), x2 = bquote(.(as.symbol(cause_var))(t))))

  ## Reshape data for plotting
  coef_res_long <- gather(coef_res_cols, key = "coef", value = "value", c_effect, c_cause)

  ## Plot coefficients
  ggplot(coef_res_long, aes(x = coef, y = value)) +
    geom_boxplot() +
    geom_violin(alpha = 0.5) +
    geom_sina(alpha = 0.5) +
    scale_x_discrete(labels = c(
      "c_effect" = c_effect_expr,
      "c_cause" = c_cause_expr
    )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    theme(text = element_text(size = 20), legend.position = "none", axis.title.x = element_blank()) +
    labs(y = "Regularized S-map Coefficient")

  ggsave(paste0(fig_date_path, "smapFig/", "coefficients_", cause_var, "_causes_", effect_var, ".tiff"),
         units = "in", width = 8, height = 8, dpi = 300, compression = 'lzw')

  ## Time series plots
  date_xaxis <- as.POSIXct(df$date, format = "%Y-%m-%d")

  # Observed and predicted plot
  ObsPred_plt <- ggplot(smap_pred_res, aes(date_xaxis)) +
    geom_line(aes(y = obs, linetype = "Observation")) +
    geom_line(aes(y = pred, linetype = "Prediction")) +
    scale_colour_grey(start = 0, end = 0.8) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    ylab("Incidence") +
    labs(linetype = effect_var) +
    theme(text = element_text(size = 20), axis.title.x = element_blank(),
          legend.position = c(0.80, 0.9), axis.title.y = element_text(angle = 0, vjust = 0.5))

  ggsave(paste0(fig_date_path, "smapFig/", "ObsPred_", cause_var, "_causes_", effect_var, ".tiff"),
         units = "in", width = 24, height = 8, dpi = 300, compression = 'lzw')

  # Coefficient plots over time
  c_effect_plot <- ggplot(coef_res_cols, aes(date_xaxis)) +
    geom_line(aes(y = c_effect)) +
    geom_hline(yintercept = quantile(na.omit(coef_res_cols$c_effect), c(0.25, 0.5, 0.75)), linetype = "dashed", color = c("gray", "black", "gray")) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    labs(y = c_effect_expr) +
    theme(text = element_text(size = 20), axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5))

  ggsave(paste0(fig_date_path, "smapFig/", "effect_", cause_var, "_causes_", effect_var, ".tiff"),
         units = "in", width = 24, height = 8, dpi = 300, compression = 'lzw')

  c_cause_plot <- ggplot(coef_res_cols, aes(date_xaxis)) +
    geom_line(aes(y = c_cause)) +
    geom_hline(yintercept = quantile(na.omit(coef_res_cols$c_cause), c(0.25, 0.5, 0.75)), linetype = "dashed", color = c("gray", "black", "gray")) +
    scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
    labs(y = c_cause_expr) +
    theme(text = element_text(size = 20), axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5))

  ggsave(paste0(fig_date_path, "smapFig/", "cause_", cause_var, "_causes_", effect_var, ".tiff"),
         units = "in", width = 24, height = 8, dpi = 300, compression = 'lzw')

  # Combine plots if needed
  merged_plot <- plot_grid(ObsPred_plt, c_effect_plot, c_cause_plot, align = "v", ncol = 1)
  merged_plot
  ggsave(paste0(fig_date_path, "smapFig/", "merged_", cause_var, "_causes_", effect_var, ".tiff"),
         units = "in", width = 24, height = 24, dpi = 300, compression = 'lzw')
}




## Function
## Use rUIC::simplex() to calculate the optimal E, output BestE_uic_simplex, detailed info in BestEforEachTP_res
UIC_simplex_res <- function(libVar, E_range, tp_range) {
  uic_simplex <- rUIC::simplex(df, lib_var = libVar, E = E_range, tau = 1, tp = tp_range, alpha = 0.05)
  BestEforEachTP_res <- list()
  for (tp_var in tp_range) {
    df_EachTP <- uic_simplex %>% dplyr::filter(tp == tp_var)
    BestEforEachTP <- with(df_EachTP, max(c(0, E[pval < 0.05])))
    BestEforEachTP_res[[tp_var + 13]] <- data.frame(df_EachTP %>% dplyr::filter(E == BestEforEachTP))
  }
  BestEforEachTP_res <- do.call(rbind, BestEforEachTP_res)
  # write.csv(BestEforEachTP_res, paste0(uic_surr_tbl_path, "result/simplexE_res_", tarVar, "_cause_", libVar, ".csv"))
  BestE_uic_simplex <- with(BestEforEachTP_res, max(c(0, E)))
  # write.csv(BestE_uic_simplex, paste0(uic_surr_tbl_path, "result/simplex_E_", tarVar, "_cause_", libVar, ".csv"))
  print(paste("E - simplex: ", libVar, BestE_uic_simplex))
}

UIC_multiSimplex_res <- function(libVar, tarVar, E_range, tp_range) {
  uic_multiSimplex <- rUIC::simplex(df, lib_var = libVar, cond_var = tarVar, E = E_range, tau = 1, tp = tp_range, alpha = 0.05)
  BestEforEachTP_multi_res <- list()
  for (tp_var in tp_range) {
    dfMulti_EachTP <- uic_multiSimplex %>% dplyr::filter(tp == tp_var)
    BestEforEachTP_multi <- with(dfMulti_EachTP, max(c(0, E[pval < 0.05])))
    BestEforEachTP_multi_res[[tp_var + 13]] <- data.frame(dfMulti_EachTP %>% dplyr::filter(E == BestEforEachTP_multi))
  }
  BestEforEachTP_multi_res <- do.call(rbind, BestEforEachTP_multi_res)
  # write.csv(BestEforEachTP_multi_res, paste0(uic_surr_tbl_path, "result/multi_simplexE_res_", tarVar, "_cause_", libVar, ".csv"))
  BestE_uic_multiSimplex <- with(BestEforEachTP_multi_res, max(c(0, E)))
  # write.csv(BestE_uic_multiSimplex, paste0(uic_surr_tbl_path, "result/multi_simplex_E_", tarVar, "_cause_", libVar, ".csv"))
  print(paste("E - multi simplex: ", libVar, BestE_uic_multiSimplex))
}


# Surr_UIC <- function(libVar, tarVar, numSurr) {
#   tp_range <- c(-12:4)
#   E_range <- c(0:20)
#   ## tarVar cause libVar?
#   uic_res <- uic.optimal(df, lib_var = libVar, tar_var = tarVar, E = E_range, tau = 1, tp = tp_range, alpha = 0.05)
#   print(paste("E - uic optimal: ", libVar, uic_res$E[1] + 1))
#
#   # UIC_simplex_res(libVar, E_range, tp_range)
#   # UIC_multiSimplex_res(libVar, tarVar, E_range, tp_range)
#
#   # compute UIC using optimal embedding dimension
#
#   write.csv(uic_res, paste0(uic_surr_tbl_path, "result/", tarVar, "_cause_", libVar, "_UIC_result.csv"))
#
#   ggplot(uic_res, aes(x = tp, y = ete)) +
#     geom_line() +
#     labs(title = paste0("UIC (", tarVar, " cause ", libVar, "?)"), shape = "pval") +
#     theme(title = element_text(face = "bold"), text = element_text(size = 20)) +
#     geom_point(aes(color = pval < 0.05), size = 4) +
#     theme(legend.position = "none") +
#     scale_color_manual(values = c("black", "red")) +
#     labs(x = "tp", y = "Effective TE", color = "p < 0.05")
#   ggsave(paste0(uic_plot_path, "uic/", tarVar, "_cause_", libVar, "_UIC.tiff"),
#          units = "in",
#          width = 10, height = 8, dpi = 300, compression = 'lzw')
#
#   ## Generate surrogate ts
#   lib_surr <- rEDM::make_surrogate_seasonal(df[[libVar]], num_surr = numSurr)
#   write.csv(lib_surr, paste0(uic_surr_tbl_path, "surrogate_dt/", tarVar, "_cause_", libVar, "_surrogate_data.csv"))
#
#   ## UIC for surrogate data
#   ### Result object
#   ete_surr <- data.frame(tp = tp_range)
#   for (i in 1:ncol(lib_surr)) {
#     block_tmp <- data.frame(effect = lib_surr[, i], cause = df[[tarVar]])
#     ete_surr_i <- uic.optimal(block_tmp,
#                               lib_var = "effect",
#                               tar_var = "cause", E = E_range, tau = 1, tp = tp_range, num_surr = 1)
#     res_ete_tp <- data.frame(ete_surr_i %>% select(tp, ete))
#     colnames(res_ete_tp) <- c("tp", paste0("ete_", i))
#     ete_surr <- merge(ete_surr, res_ete_tp, by = "tp")
#     # tp & ete: [tp's range] x [len(surr) +1]
#   }
#   write.csv(ete_surr, paste0(uic_surr_tbl_path, "result/", tarVar, "_cause_", libVar, "_ete_tp_surrogate_data.csv"))
#
#   ## Calculate 95% CI
#   quantile_list <- list()
#   for (i in 1:nrow(ete_surr)) {
#     ete_surr_noTP <- subset(ete_surr, select = -c(tp))
#     res_quantile <- quantile(unlist(ete_surr_noTP[i,]), probs = c(0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99))
#     quantile_list[[paste0("tp_", i - 9)]] <- res_quantile
#   }
#   df_quantiles <- data.frame(do.call(rbind, quantile_list))
#   write.csv(df_quantiles, paste0(uic_surr_tbl_path, "result/", tarVar, "_cause_", libVar, "_quantile_surrogate_data.csv"))
#
#   significant_UIC <- data.frame(tp = uic_res$tp,
#                                 ete = uic_res$ete,
#                                 pval = uic_res$pval,
#                                 quantile_90 = df_quantiles$X90.,
#                                 quantile_95 = df_quantiles$X95.,
#                                 quantile_97.5 = df_quantiles$X97.5.,
#                                 quantile_99 = df_quantiles$X99.)
#   write.csv(significant_UIC, paste0(uic_surr_tbl_path, tarVar, "_cause_", libVar, "_uic_p_surr_res.csv"))
#
#   ggplot(significant_UIC, aes(x = tp)) +
#     geom_line(aes(y = ete), linetype = "solid") +
#     geom_line(aes(y = quantile_95), linetype = "longdash") +
#     geom_line(aes(y = quantile_90), linetype = "dotted") +
#     labs(title = paste0("UIC (", tarVar, " cause ", libVar, "?)"), shape = "pval") +
#     theme(title = element_text(face = "bold"), text = element_text(size = 20)) +
#     geom_point(aes(y = ete, color = pval < 0.05), size = 4) +
#     theme(legend.position = "none") +
#     scale_color_manual(values = c("black", "red")) +
#     labs(x = "tp", y = "Effective TE", color = "p < 0.05")
#
#   ggsave(paste0(uic_plot_path, "uic_surr/", tarVar, "_cause_", libVar, "_UIC_surr.tiff"),
#          units = "in",
#          width = 10, height = 8, dpi = 300, compression = 'lzw')
# }

# # use saving path from upper uic save path
# reg_smap_func <- function(lib_var, tar_var) {
#   print(paste("effected:",lib_var,"cause",tar_var))
#   ## calculate the Best E
#   tp_range <- c(-12:4)
#   E_range <- c(0:20)
#   ## tarVar cause libVar?
#   # uic_res <- uic.optimal(df, lib_var = lib_var, tar_var = tar_var, E = E_range, tau = 1, tp = tp_range)
#   # BestE <- uic_res$E[1]
#
#   ## UIC univariate simplex
#   # uic_simplex <- rUIC::simplex(df, lib_var = lib_var, E = E_range, tau = 1, tp = tp_range, alpha = 0.05)
#   # BestEforEachTP_res <- list()
#   # for (tp_var in tp_range) {
#   #   df_EachTP <- uic_simplex %>% dplyr::filter(tp == tp_var)
#   #   BestEforEachTP <- with(df_EachTP, max(c(0, E[pval < 0.05])))
#   #   BestEforEachTP_res[[tp_var + 13]] <- data.frame(df_EachTP %>% dplyr::filter(E == BestEforEachTP))
#   # }
#   # BestEforEachTP_res <- do.call(rbind, BestEforEachTP_res)
#   # BestE <- with(BestEforEachTP_res, max(c(0, E)))
#   # print(paste("E - uic simplex: ", lib_var, BestE))
#
#   ## tarVar cause libVar?
#   uic_res <- uic.optimal(df, lib_var = lib_var, tar_var = tar_var, E = E_range, tau = 1, tp = tp_range, alpha = 0.05)
#   BestE <- uic_res$E[1] + 1
#   print(paste("E - uic optimal: ", lib_var, BestE))
#
#   lag_block_target <- make_block(df[[tar_var]], max_lag = BestE)
#
#   # if have multiple lib_var
#   for (var in 1:length(lib_var)) {
#     uic_res <- read.csv(paste0(uic_surr_tbl_path, tar_var, "_cause_", lib_var[[var]], "_uic_p_surr_res.csv"))
#     uic_res_pval <- uic_res[uic_res$pval < 0.05,]
#     BestTP <- uic_res_pval[which.max(uic_res_pval$ete),]$tp
#     print(paste("Best tp:", BestTP))
#
#     if (BestTP == 0) {
#       lag_block_lib <- dplyr::lead(df[[lib_var]], n = 1) # tp = 0
#     }else {
#       lag_block_lib <- dplyr::lag(df[[lib_var]], n = abs(BestTP) - 1)
#     }
#   }
#
#   smap_block <- cbind(lag_block_target,
#                       lag_block_lib)
#   smap_block <- smap_block[, 2:ncol(smap_block)]
#   write.csv(smap_block, paste0(res_save_path, "table/smap/block/", lib_var, "(effected)", tar_var, "(cause_var)"))
#
#   theta_range <- c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
#   stat_res <- list()
#   for (theta in theta_range) {
#     print(theta)
#     rsmap_ridge <- macamts::extended_lnlp(smap_block, theta = theta, lambda = 0.1,
#                                           regularized = T, alpha = 0, random_seed = 1234, save_smap_coefficients = TRUE)
#     stat_res <- rbind(stat_res, rsmap_ridge$stats)
#   }
#   stat_res <- cbind(theta_range, stat_res)
#   BestTheta <- stat_res[which.min(stat_res$rmse), "theta_range"]
#   print(paste("theta - reg. s-map: ", lib_var, tar_var, BestTheta))
#
#   rsmap_ridge <- macamts::extended_lnlp(smap_block, theta = BestTheta, lambda = 0.1,
#                                         regularized = T, alpha = 0, random_seed = 1234, save_smap_coefficients = TRUE)
#   smap_pred_res <- rsmap_ridge$model_output
#   write.csv(smap_pred_res, paste0(res_save_path, "table/smap/pred_res/", lib_var, "(effected)", tar_var, "(cause_var)"))
#   coef_res <- rsmap_ridge$smap_coefficients
#   write.csv(coef_res, paste0(res_save_path, "table/smap/coef/", lib_var, "(effected)", tar_var, "(cause_var)"))
#
#   maxValue <- max(max(na.omit(smap_pred_res$obs)), max(na.omit(smap_pred_res$pred))) + 1
#   ggplot(smap_pred_res, aes(x = obs, y = pred, color = time)) +
#     geom_abline(slope = 1, linetype = "dashed", color = "black") +
#     geom_point() +
#     # same xlim, ylim
#     xlim(c(NA, maxValue)) +
#     ylim(c(NA, maxValue))
#   ggsave(paste0(fig_date_path, "smapFig/", "pred_obs_",tar_var, "(cause_var)", lib_var, "(effected)", ".tiff"),
#          units = "in",
#          width = 8, height = 8, dpi = 300, compression = 'lzw')
#
#   coef_res_cols <- data.frame(time = coef_res$time,
#                               c_1 = coef_res$c_1,
#                               c_2 = coef_res[[paste0("c_", BestE + 1)]])
#   write.csv(coef_res_cols, paste0(res_save_path, "table/smap/coef/","sumRes_", lib_var, "(effected)", tar_var, "(cause_var)"))
#
#   c1_var <- tar_var
#   c1_math_expression <- substitute(frac(partialdiff * x1, partialdiff * x2), list(x1 = bquote(.(as.symbol(c1_var))(t + 1)), x2 = bquote(.(as.symbol(c1_var))(t))))
#   c2_var <- lib_var
#   c2_math_expression <- substitute(frac(partialdiff * x1, partialdiff * x2), list(x1 = bquote(.(as.symbol(c1_var))(t + 1)), x2 = bquote(.(as.symbol(c2_var))(t))))
#
#
#   coef_res_long <- gather(coef_res, key = "coef", value = "value", c_1, c_2)
#
#   # sinaplot, boxplot
#   ggplot(coef_res_long, aes(x = coef, y = value)) +
#     geom_boxplot() +
#     geom_violin(alpha = .5) +
#     geom_sina(alpha = .5) +
#
#     scale_x_discrete(labels = c(
#       "c_1" = c1_math_expression,
#       "c_2" = c2_math_expression
#     )) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
#     ## mean
#     geom_hline(yintercept = mean(na.omit(coef_res$c_1)), linetype = "dashed", color = "gray") +
#     geom_hline(yintercept = mean(na.omit(coef_res$c_2)), linetype = "dashed", color = "gray") +
#
#     theme(text = element_text(size = 20), legend.position = c(.25, .65), axis.title.x = element_blank()) +
#     labs(x = paste0("Interactions Between ", c1_var, " and ", c2_var, y = "S-map coefficients"), y = "Regularized S-map Coefficient")
#   ggsave(paste0(fig_date_path, "smapFig/","boxplot_" , "pred_obs_",tar_var, "(cause_var)", lib_var, "(effected)", ".tiff"),
#          units = "in",
#          width = 8, height = 8, dpi = 300, compression = 'lzw')
#
#   date_xaxis <- as.POSIXct(df$Date, format = "%Y-%m-%d")
#
#   # ts of observed and predicted
#   ObsPred_plt <- ggplot(smap_pred_res, aes(date_xaxis)) +
#     geom_line(aes(y = obs, linetype = "Observation")) +
#     geom_line(aes(y = pred, linetype = "Prediction")) +
#     scale_colour_grey(start = 0, end = .8) +
#     scale_x_datetime(
#       date_breaks = "1 year",
#       date_labels = "%Y",
#       minor_breaks = "1 month"
#     ) +
#     ylab("Incidence") +
#     labs(linetype = c1_var) +
#     theme(text = element_text(size = 20), axis.title.x = element_blank(), legend.position = c(.80, .9),
#           axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
#   ObsPred_plt
#   ggsave(paste0(fig_date_path, "smapFig/","ObsPred_plt_" , "pred_obs_",tar_var, "(cause_var)", lib_var, "(effected)", ".tiff"),
#          units = "in",
#          width = 24, height = 8, dpi = 300, compression = 'lzw')
#
#   # plot of c_1
#   c1Y_plt <- ggplot(coef_res, aes(date_xaxis)) +
#     geom_hline(yintercept = quantile(na.omit(coef_res$c_1), 0.25), linetype = "dashed", color = "gray") +
#     geom_hline(yintercept = quantile(na.omit(coef_res$c_1), 0.5), linetype = "dashed", color = "black") +
#     geom_hline(yintercept = quantile(na.omit(coef_res$c_1), 0.75), linetype = "dashed", color = "gray") +
#     geom_line(aes(y = c_1)) +
#     # geom_line(aes(y = c_2, linetype = "c_2")) +
#     scale_colour_grey(start = 0, end = .8) +
#     scale_x_datetime(
#       date_breaks = "1 year",
#       date_labels = "%Y",
#       minor_breaks = "1 month"
#     ) +
#     labs(y = c1_math_expression) +
#     theme(text = element_text(size = 20), axis.title.x = element_blank(),
#           axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
#   c1Y_plt
#   ggsave(paste0(fig_date_path, "smapFig/","c1Y_plt_" , "pred_obs_",tar_var, "(cause_var)", lib_var, "(effected)", ".tiff"),
#          units = "in",
#          width = 24, height = 8, dpi = 300, compression = 'lzw')
#
#   # plot of c_2
#   c2X_plt <- ggplot(coef_res, aes(date_xaxis)) +
#     geom_hline(yintercept = quantile(na.omit(coef_res$c_2), 0.25), linetype = "dashed", color = "gray") +
#     geom_hline(yintercept = quantile(na.omit(coef_res$c_2), 0.5), linetype = "dashed", color = "black") +
#     geom_hline(yintercept = quantile(na.omit(coef_res$c_2), 0.75), linetype = "dashed", color = "gray") +
#     geom_line(aes(y = c_2)) +
#     scale_colour_grey(start = 0, end = .8) +
#     labs(y = c2_math_expression) +
#     scale_x_datetime(
#       date_breaks = "1 year",
#       date_labels = "%Y",
#       minor_breaks = "1 month"
#     ) +
#     # geom_hline(yintercept = cor_result, linetype = "dashed", color = "red") +
#     theme(text = element_text(size = 20), axis.title.x = element_blank(),
#           axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
#   c2X_plt
#   ggsave(paste0(fig_date_path, "smapFig/","c2X_plt_" , "pred_obs_",tar_var, "(cause_var)", lib_var, "(effected)", ".tiff"),
#          units = "in",
#          width = 24, height = 8, dpi = 300, compression = 'lzw')
#
#   # merge 3 plots
#   merged_plt <- plot_grid(ObsPred_plt, c1Y_plt, c2X_plt, align = "v", ncol = 1)
#   merged_plt
#   ggsave(paste0(fig_date_path, "smapFig/","merged_plt_" , "pred_obs_",tar_var, "(cause_var)", lib_var, "(effected)", ".tiff"),
#          units = "in",
#          width = 24, height = 10, dpi = 300, compression = 'lzw')
# }