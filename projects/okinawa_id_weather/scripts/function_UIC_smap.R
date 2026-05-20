run_uic_smap_pair <- function(df_log,
                              effect_var,
                              cause_var,
                              uic_surr_tbl_path,
                              res_save_path,
                              xmap_plot_path,
                              fig_date_path) {
print(paste("Effected variable:", effect_var, "| Cause variable:", cause_var))
df <- df_log

## Calculate the optimal embedding dimension (E)
tp_range <- -12:0
E_range <- 0:20
uic_res <- uic.optimal(df, lib_var = effect_var, tar_var = cause_var, E = E_range, tau = 1, tp = tp_range, alpha = 0.05)
BestE <- uic_res$E[1] + 1
print(paste("Optimal E for", effect_var, ":", BestE))

## Prepare lagged variables based on the best time lag (tp)
# uic_surr_tbl_path <- "result/FluSub_JP/2024-10-08/table/uic/"
uic_res <- read.csv(paste0(uic_surr_tbl_path, cause_var, "_cause_", effect_var, "_uic_p_surr_res.csv"))
uic_res_pval <- uic_res[uic_res$pval < 0.05,]
BestTP <- uic_res_pval[which.max(uic_res_pval$ete),]$tp
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
date_xaxis <- as.POSIXct(df$Date, format = "%Y-%m-%d")

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
