# flu_subtype_mdr_smap_primary_pipeline.R
# Primary EDM pipeline for influenza subtype interaction analysis using UIC + MDR S-map.
#
# Purpose:
#   1) Detect directional, lagged subtype interactions by UIC with seasonal surrogate correction.
#   2) Use surrogate-significant UIC links to construct multiview-distance blocks.
#   3) Estimate interaction polarity/strength with MDR S-map as the primary S-map analysis.
#   4) Export manuscript-ready tables, Figure 3/4/5 inputs, and figure files.
#
# Expected input:
#   data/FluSub_jp/FluSub_11to19_jp_per_20240925.csv
#   Required columns: year, week, B, A_H1N1, A_H3N2
#
# Notes:
#   - Run from the project root where data/ and result/ exist.
#   - For a quick dry run, set config$num_surr <- 50. For final results use 2000.
#   - MDR functions are taken from either macamts or macam, whichever is installed.

rm(list = ls())

# -----------------------------
# 0. Bootstrap and packages
# -----------------------------
setwd("~/Desktop/_GSAIS_/mzmtlab/microbiome dynamics/edm_code/projects/jp_flu_subtypes/scripts")
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
setwd(workspace_root())
rm(source_edm_bootstrap)

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, # データ加工
  tidyr, # データの縦横変換
  purrr, # 関数型処理、反復処理
  tibble, # データフレームの拡張
  readr, # csvなどの読み込み
  stringr, # 文字列処理
  forcats, # カテゴリ変数の処理
  scales, # ggplot の軸・ラベル調整
  lubridate, # 日付・時刻の操作
  ISOweek, # ISO週番号の計算
  ggplot2, # グラフ作成
  cowplot, # 複数図の結合
  sinaplot, # データの分布プロット（バイオリンプロット＋点プロット）
  rUIC, # Unified Information Criterion 関連
  rEDM, # Empirical Dynamic Modeling
  grid # 図の細かい制御
)
# Load library
packageVersion("rEDM") # v 0.7.5
packageVersion("macamts") # v 0.1.4
packageVersion("rUIC") # v 0.9.12

resolve_mdr_pkg <- function() {
  candidates <- c("macamts", "macam")
  for (pkg in candidates) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      return(pkg)
    }
  }
  stop(
    "Neither 'macamts' nor 'macam' is installed. Install one of them, e.g. remotes::install_github('ong8181/macam').",
    call. = FALSE
  )
}

mdr_pkg <- resolve_mdr_pkg()
message("Using MDR package: ", mdr_pkg)

mdr_fun <- function(fun) {
  getExportedValue(mdr_pkg, fun)
}

# Compatibility wrapper: macam/macamts >= 0.2.2 uses max_lag; older macamts uses E_effect_var.
# バージョン差を吸収する wrapper
make_block_mvd_compat <- function(block, uic_res, effect_var, max_lag, ...) {
  f <- mdr_fun("make_block_mvd")
  dots <- list(...)
  formal_names <- names(formals(f))

  args <- c(list(block = block, uic_res = uic_res, effect_var = effect_var), dots)
  if (!is.null(dots$tp_adjust) && !("tp_adjust" %in% formal_names)) {
    warning("Installed make_block_mvd() has no tp_adjust argument. UIC lag and one-step MDR target may be offset by one week. Update macam/macamts if possible.")
  }
  if ("max_lag" %in% formal_names) {
    args$max_lag <- max_lag
  } else if ("E_effect_var" %in% formal_names) {
    args$E_effect_var <- max_lag
  } else {
    warning("make_block_mvd() has neither max_lag nor E_effect_var; using package default for effect-variable lags.")
  }

  # Drop arguments unsupported by the installed package version.
  args <- args[names(args) %in% formal_names]
  do.call(f, args)
}

call_mdr_function <- function(fun, ...) {
  f <- mdr_fun(fun) # MDR パッケージから関数を取得する
  args <- list(...) # 渡された引数をリスト化する
  formal_names <- names(formals(f)) # 関数の引数名を取得する
  args <- args[names(args) %in% formal_names] # 対応していない引数を削除する
  do.call(f, args) # 関数を実行する
}

# MVD 計算結果の形式を統一
compute_mvd_compat <- function(...) {
  out <- call_mdr_function("compute_mvd", ...)
  if (is.matrix(out) || inherits(out, "dist")) {
    return(list(multiview_dist = as.matrix(out)))
  }
  if (is.list(out) && "multiview_dist" %in% names(out)) {
    return(out)
  }
  stop("compute_mvd() returned an unsupported object. Inspect the installed MDR package version.", call. = FALSE)
}

s_map_mdr_compat <- function(...) {
  call_mdr_function("s_map_mdr", ...)
}

# -----------------------------
# 1. Configuration
# -----------------------------
config <- list(
  data_file = "data/FluSub_jp/FluSub_11to19_jp_per_20240925.csv",
  out_dir = file.path("result", "FluSub_JP", paste0(format(Sys.Date(), "%Y%m%d"), "_UIC_MDR_primary")),
  subtype_vars = c("B", "A_H1N1", "A_H3N2"),

  # UIC settings. Primary analysis excludes tp = 0 to avoid contemporaneous seasonal synchrony.
  E_range = 1:20, # 埋め込み次元またはラグ数候補の範囲(過去何ステップ分の情報を状態空間再構成に使うか、または対象変数の自己履歴をどの程度含めるかに関係します)
  tp_range = -12:-1, # 予測ラグの範囲(例えば、 tp=-3 なら「A/H1N1 の 3週前の値」が、現在または1週後の B に関係するか)
  tau = 1, # 予測ステップ数（デフォルトは1, 1週間刻みでラグを取る）
  alpha = 0.05, # 有意水準
  num_surr = 2000, # サロゲートデータの数
  season_period = 52, # 季節周期（季節周期を 52 週に設定、単位はステップ数、日本でのインフルエンザシーズンを考慮）
  random_seed = 1234, # 乱数シード

  # MDR S-map settings.
  smap_tp = 1, # S-map の予測ターゲットとなる時間ステップ（1週間後の値を予測）
  mdr_include_var = "strongest_only", # 各効果変数に対して、UICで最も強い関係を示した原因変数だけを MDR S-map に入れる設定
  mdr_E = 3, # MDR S-map に使用する埋め込み次元：状態空間を構成する際に3次元の情報を使う
  n_ssr = 2000, # サロゲートデータの数
  k = NULL, # if NULL, floor(sqrt(n_ssr)) is used. 近傍点数？
  theta_grid = c(
    0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2,
    0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8
  ), # 拡散係数のグリッド値、S-mapでは、theta が大きいほど、近い状態にある点を重視する局所モデルになる

  # 正則化付き MDR S-map を使うかどうか。感度分析として、TRUE のままで解析を進めることも可能
  ridge_regularized_mdr = FALSE, # 主解析は使わない
  lambda_grid = c(0), # 正則化の強さのグリッド値。0は正則化なし
  alpha_glmnet = 0, # alpha = 0：Ridge, alpha = 1：Lasso, 0 < alpha < 1：Elastic Net

  # Figure export.
  dpi = 300,
  fig_width = 12,
  fig_height = 8
)
if (is.null(config$k)) config$k <- floor(sqrt(config$n_ssr))

# -----------------------------
# 2. Output paths
# -----------------------------
safe_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

paths <- list(
  root = safe_dir(config$out_dir),
  tables = safe_dir(file.path(config$out_dir, "tables")),
  uic_tables = safe_dir(file.path(config$out_dir, "tables", "uic")),
  mdr_tables = safe_dir(file.path(config$out_dir, "tables", "mdr_smap")),
  fig = safe_dir(file.path(config$out_dir, "figures")),
  uic_fig = safe_dir(file.path(config$out_dir, "figures", "uic")),
  mdr_fig = safe_dir(file.path(config$out_dir, "figures", "mdr_smap")),
  manuscript = safe_dir(file.path(config$out_dir, "manuscript_ready"))
)

save_csv <- function(x, path) {
  readr::write_csv(x, path, na = "")
  invisible(path)
}

save_rds <- function(x, path) {
  saveRDS(x, path)
  invisible(path)
}

subtype_label <- function(x) {
  dplyr::recode(x,
    "A_H1N1" = "A/H1N1",
    "A_H3N2" = "A/H3N2",
    "B" = "Type B",
    .default = x
  )
}

# 因果方向・影響方向のラベルを作る関数
edge_label <- function(cause, effect, lag_weeks = NULL) {
  base <- paste0(subtype_label(cause), " \u2192 ", subtype_label(effect))
  if (is.null(lag_weeks)) {
    return(base)
  }
  paste0(base, " (", lag_weeks, " wk)")
}

theme_pub <- function(base_size = 14) {
  cowplot::theme_cowplot(font_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}

theme_set(theme_pub())

# -----------------------------
# 3. Data preparation
# -----------------------------
read_prepare_flu <- function(data_file, subtype_vars) {
  if (!file.exists(data_file)) {
    stop("Data file not found: ", data_file, call. = FALSE)
  }

  raw <- read.csv(data_file, check.names = FALSE)

  required <- c("year", "week", subtype_vars)
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  raw$Date <- ISOweek::ISOweek2date(
    paste0(
      raw$year,
      "-W",
      sprintf("%02d", as.integer(raw$week)),
      "-1"
    )
  )

  raw$Date <- as.Date(raw$Date)
  raw <- raw[order(raw$Date), ]

  out <- data.frame(Date = raw$Date)

  for (v in subtype_vars) {
    x <- suppressWarnings(as.numeric(raw[[v]]))

    if (anyNA(x)) {
      stop(
        "Subtype variable contains non-numeric or missing values: ",
        v,
        call. = FALSE
      )
    }

    out[[v]] <- log(x + 1)
  }

  keep <- c(
    TRUE,
    vapply(out[subtype_vars], function(x) sum(x, na.rm = TRUE) > 0, logical(1))
  )

  out <- out[, keep, drop = FALSE]

  if (anyNA(out)) {
    stop(
      "Prepared data contains NA values. Please inspect missing weeks or subtype values.",
      call. = FALSE
    )
  }

  out
}

df_log <- read_prepare_flu(config$data_file, config$subtype_vars)
df_model <- df_log %>%
  dplyr::select(-Date) %>%
  as.data.frame()
vars <- names(df_model)

save_csv(df_log, file.path(paths$tables, "prepared_log1p_timeseries.csv"))

# -----------------------------
# 4. Seasonal surrogate for weekly data
# -----------------------------
make_seasonal_surrogates <- function(ts, num_surr, period = 52, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ts <- as.numeric(ts)
  if (any(!is.finite(ts))) stop("Input time series contains non-finite values.", call. = FALSE)

  n <- length(ts)
  season_index <- rep(seq_len(period), length.out = n)

  seasonal_fit <- smooth.spline(
    x = c(season_index - period, season_index, season_index + period),
    y = c(ts, ts, ts)
  )
  seasonal_cycle <- predict(seasonal_fit, season_index)$y
  residuals <- ts - seasonal_cycle

  out <- matrix(NA_real_, nrow = n, ncol = num_surr)
  for (i in seq_len(num_surr)) {
    out[, i] <- seasonal_cycle + sample(residuals, size = n, replace = FALSE)
  }
  colnames(out) <- paste0("surr_", seq_len(num_surr))
  out
}

# -----------------------------
# 5. UIC with seasonal surrogate correction
# -----------------------------
run_uic_pair <- function(df_model, effect_var, cause_var, config, paths) {
  message("UIC: ", cause_var, " -> ", effect_var)

  obs <- rUIC::uic.optimal(
    as.data.frame(df_model),
    lib_var = effect_var,
    tar_var = cause_var,
    E = config$E_range,
    tau = config$tau,
    tp = config$tp_range,
    alpha = config$alpha
  ) %>%
    as_tibble() %>%
    filter(tp %in% config$tp_range) %>%
    arrange(tp)

  if (!("ete" %in% names(obs))) {
    stop("rUIC::uic.optimal() output does not include 'ete'.", call. = FALSE)
  }
  if (!("te" %in% names(obs))) {
    obs <- obs %>% mutate(te = .data$ete)
  }
  if (!("pval" %in% names(obs))) {
    obs <- obs %>% mutate(pval = NA_real_)
  }

  effect_surr <- make_seasonal_surrogates(
    df_model[[effect_var]],
    num_surr = config$num_surr,
    period = config$season_period,
    seed = config$random_seed
  )

  surr_ete <- matrix(NA_real_, nrow = nrow(obs), ncol = config$num_surr)
  rownames(surr_ete) <- paste0("tp_", obs$tp)
  colnames(surr_ete) <- paste0("surr_", seq_len(config$num_surr))

  for (i in seq_len(config$num_surr)) {
    tmp <- data.frame(effect = effect_surr[, i], cause = df_model[[cause_var]])
    sres <- rUIC::uic.optimal(
      as.data.frame(tmp),
      lib_var = "effect",
      tar_var = "cause",
      E = config$E_range,
      tau = config$tau,
      tp = config$tp_range,
      alpha = config$alpha
    ) %>%
      as_tibble() %>%
      select(tp, ete)
    surr_ete[, i] <- sres$ete[match(obs$tp, sres$tp)]
  }

  qfun <- function(x, p) as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
  q90 <- apply(surr_ete, 1, qfun, p = 0.90)
  q95 <- apply(surr_ete, 1, qfun, p = 0.95)
  q975 <- apply(surr_ete, 1, qfun, p = 0.975)
  q99 <- apply(surr_ete, 1, qfun, p = 0.99)

  p_emp <- vapply(seq_len(nrow(obs)), function(r) {
    (1 + sum(surr_ete[r, ] >= obs$ete[r], na.rm = TRUE)) / (1 + config$num_surr)
  }, numeric(1))

  max_surr <- apply(surr_ete, 2, max, na.rm = TRUE)
  q95_global <- qfun(max_surr, 0.95)
  p_global <- vapply(obs$ete, function(x) {
    (1 + sum(max_surr >= x, na.rm = TRUE)) / (1 + config$num_surr)
  }, numeric(1))

  res <- obs %>%
    mutate(
      effect_var = effect_var,
      cause_var = cause_var,
      q90 = q90,
      q95 = q95,
      q975 = q975,
      q99 = q99,
      q95_global = q95_global,
      p_emp = p_emp,
      p_global = p_global,
      sig_pointwise = .data$p_emp < config$alpha & .data$ete > .data$q95,
      sig_global = .data$p_global < config$alpha & .data$ete > .data$q95_global
    )

  tag <- paste0(cause_var, "_to_", effect_var)
  save_csv(res, file.path(paths$uic_tables, paste0(tag, "_uic_surrogate_corrected.csv")))
  save_csv(as.data.frame(surr_ete), file.path(paths$uic_tables, paste0(tag, "_surrogate_ete_matrix.csv")))

  p <- ggplot(res, aes(x = .data$tp)) +
    geom_line(aes(y = .data$ete), linewidth = 0.6) +
    geom_line(aes(y = .data$q95), linetype = "longdash", linewidth = 0.5) +
    geom_line(aes(y = .data$q90), linetype = "dotted", linewidth = 0.5) +
    geom_hline(yintercept = q95_global, linetype = "dotdash", linewidth = 0.5) +
    geom_point(aes(y = .data$ete, shape = .data$sig_global), size = 2.8) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19)) +
    scale_x_continuous(breaks = config$tp_range) +
    labs(
      title = edge_label(cause_var, effect_var),
      x = "UIC tp; negative values mean that the cause precedes the effect",
      y = "Effective transfer entropy",
      shape = "surrogate p < 0.05"
    ) +
    theme_pub(13)

  ggsave(file.path(paths$uic_fig, paste0(tag, "_Figure3_component_UIC.tiff")), p,
    width = 8, height = 6, dpi = config$dpi, compression = "lzw"
  )

  res
}

uic_all <- purrr::map_dfr(vars, function(effect_var) {
  purrr::map_dfr(setdiff(vars, effect_var), function(cause_var) {
    run_uic_pair(df_model, effect_var, cause_var, config, paths)
  })
})

uic_all <- uic_all %>%
  mutate(
    p_global_fdr = p.adjust(.data$p_global, method = "BH"),
    p_emp_fdr = p.adjust(.data$p_emp, method = "BH"),
    sig_primary = .data$p_global_fdr < config$alpha & .data$ete > .data$q95_global
  )

save_csv(uic_all, file.path(paths$tables, "uic_all_pairs_surrogate_corrected.csv"))

selected_links <- uic_all %>%
  filter(.data$sig_primary, .data$tp < 0) %>%
  group_by(.data$effect_var, .data$cause_var) %>%
  arrange(desc(.data$ete), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    lag_weeks = abs(.data$tp),
    edge = edge_label(.data$cause_var, .data$effect_var, .data$lag_weeks),
    pval_raw_uic = .data$pval,
    pval = .data$p_global_fdr,
    te = ifelse(is.na(.data$te), .data$ete, .data$te)
  )

save_csv(selected_links, file.path(paths$tables, "uic_selected_links_for_MDR.csv"))

if (nrow(selected_links) == 0) {
  stop("No UIC links passed the surrogate/FDR primary criterion. Inspect uic_all_pairs_surrogate_corrected.csv.", call. = FALSE)
}

# Combined UIC figure for manuscript Figure 3.
uic_plot_data <- uic_all %>%
  semi_join(selected_links %>% select(effect_var, cause_var), by = c("effect_var", "cause_var")) %>%
  mutate(edge = edge_label(.data$cause_var, .data$effect_var))

fig3 <- ggplot(uic_plot_data, aes(x = .data$tp)) +
  geom_line(aes(y = .data$ete), linewidth = 0.6) +
  geom_line(aes(y = .data$q95), linetype = "longdash", linewidth = 0.4) +
  geom_line(aes(y = .data$q90), linetype = "dotted", linewidth = 0.4) +
  geom_hline(aes(yintercept = .data$q95_global), linetype = "dotdash", linewidth = 0.35) +
  geom_point(aes(y = .data$ete, shape = .data$sig_primary), size = 2) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19)) +
  scale_x_continuous(breaks = config$tp_range) +
  facet_wrap(~edge, scales = "free_y") +
  labs(
    x = "UIC tp; negative values mean that the cause precedes the effect",
    y = "Effective transfer entropy",
    shape = "primary significant",
    title = "Figure 3. UIC with 52-week seasonal surrogate thresholds"
  ) +
  theme_pub(13)

ggsave(file.path(paths$fig, "Figure3_UIC_seasonal_surrogate_selected_links.tiff"), fig3,
  width = 13, height = 7, dpi = config$dpi, compression = "lzw"
)

# -----------------------------
# 6. MDR S-map per effect variable
# -----------------------------
coef_name_map <- function(block_mvd, coef_df) {
  # In macam/macamts extended_lnlp, c_1 corresponds to the first block column,
  # c_2 to the second block column, etc. c_0, if present, is the intercept.
  tibble(
    block_col = names(block_mvd),
    coef_col = paste0("c_", seq_along(names(block_mvd))),
    exists = paste0("c_", seq_along(names(block_mvd))) %in% names(coef_df)
  )
}

coef_dates <- function(coef_df, df_log) {
  if ("time" %in% names(coef_df)) {
    idx <- suppressWarnings(as.integer(coef_df$time))
  } else {
    idx <- seq_len(nrow(coef_df))
  }
  idx[!is.finite(idx) | idx < 1 | idx > nrow(df_log)] <- NA_integer_
  df_log$Date[idx]
}

pred_dates <- function(pred_df, df_log) {
  if ("time" %in% names(pred_df)) {
    idx <- suppressWarnings(as.integer(pred_df$time))
  } else {
    idx <- seq_len(nrow(pred_df))
  }
  idx[!is.finite(idx) | idx < 1 | idx > nrow(df_log)] <- NA_integer_
  df_log$Date[idx]
}

run_mdr_for_effect <- function(effect_var, selected_links, df_model, df_log, config, paths) {
  message("MDR S-map for effect: ", effect_var)

  links_eff <- selected_links %>% filter(.data$effect_var == effect_var)
  if (nrow(links_eff) == 0) {
    return(NULL)
  }

  # Estimate effect-variable embedding dimension for the effect-history part of block_mvd.
  simp_x <- rUIC::simplex(df_model, lib_var = effect_var, E = config$E_range, tp = config$smap_tp)
  simp_x <- as_tibble(simp_x)
  Ex <- simp_x[which.min(simp_x$rmse), "E", drop = TRUE]
  Ex <- max(1, as.integer(Ex))
  save_csv(simp_x, file.path(paths$mdr_tables, paste0(effect_var, "_simplex_for_effect_embedding.csv")))

  # Critical: make_block_mvd() internally filters on the column named 'pval'.
  # Therefore pval is replaced by the surrogate/FDR p-value for MDR block construction.
  uic_for_block <- links_eff %>%
    transmute(
      effect_var = .data$effect_var,
      cause_var = .data$cause_var,
      E = .data$E,
      tp = .data$tp,
      te = .data$te,
      ete = .data$ete,
      pval_raw_uic = .data$pval_raw_uic,
      p_global = .data$p_global,
      p_global_fdr = .data$p_global_fdr,
      pval = .data$p_global_fdr,
      q95_global = .data$q95_global,
      lag_weeks = .data$lag_weeks,
      edge = .data$edge
    )

  save_csv(uic_for_block, file.path(paths$mdr_tables, paste0(effect_var, "_uic_links_entering_MDR.csv")))

  # tp_adjust keeps the UIC lag interpretation consistent with the one-step S-map.
  # Example: UIC tp = -9 and S-map tp = 1 -> use cause_tp-8, so the target effect(t+1)
  # is paired with cause(t-8), preserving a 9-week cause-to-effect lag.
  block_mvd <- make_block_mvd_compat(
    block = df_model,
    uic_res = as.data.frame(uic_for_block),
    effect_var = effect_var,
    max_lag = Ex,
    cause_var_colname = "cause_var",
    include_var = config$mdr_include_var,
    p_threshold = config$alpha,
    tp_adjust = config$smap_tp,
    sort_tp = TRUE,
    silent = FALSE
  )

  save_csv(as_tibble(block_mvd), file.path(paths$mdr_tables, paste0(effect_var, "_block_mvd.csv")))

  mvd_E <- min(config$mdr_E, ncol(block_mvd))
  if (mvd_E < 1) stop("mvd_E became < 1; inspect block_mvd.", call. = FALSE)

  dist_info <- compute_mvd_compat(
    block_mvd = block_mvd,
    effect_var = effect_var,
    E = mvd_E,
    tp = config$smap_tp,
    n_ssr = config$n_ssr,
    k = config$k,
    random_seed = config$random_seed,
    distance_only = FALSE,
    silent = FALSE
  )
  multiview_dist <- dist_info$multiview_dist
  save_rds(dist_info, file.path(paths$mdr_tables, paste0(effect_var, "_multiview_distance_and_embeddings.rds")))

  param_grid <- expand.grid(
    theta = config$theta_grid,
    lambda = config$lambda_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  param_res <- purrr::pmap_dfr(param_grid, function(theta, lambda) {
    fit <- s_map_mdr_compat(
      block_mvd = block_mvd,
      dist_w = multiview_dist,
      tp = config$smap_tp,
      theta = theta,
      regularized = config$ridge_regularized_mdr,
      lambda = lambda,
      alpha = config$alpha_glmnet,
      save_smap_coefficients = FALSE,
      random_seed = config$random_seed
    )
    as_tibble(fit$stats) %>% mutate(theta = theta, lambda = lambda)
  }) %>%
    select(any_of(c("N", "theta", "lambda", "rho", "mae", "rmse")), everything())

  save_csv(param_res, file.path(paths$mdr_tables, paste0(effect_var, "_MDR_parameter_search.csv")))

  best <- param_res %>%
    filter(is.finite(.data$rmse)) %>%
    arrange(.data$rmse, desc(.data$rho)) %>%
    slice(1)
  if (nrow(best) == 0) stop("No valid MDR parameter result for ", effect_var, call. = FALSE)

  final_fit <- s_map_mdr_compat(
    block_mvd = block_mvd,
    dist_w = multiview_dist,
    tp = config$smap_tp,
    theta = best$theta[[1]],
    regularized = config$ridge_regularized_mdr,
    lambda = best$lambda[[1]],
    alpha = config$alpha_glmnet,
    save_smap_coefficients = TRUE,
    random_seed = config$random_seed
  )

  save_rds(final_fit, file.path(paths$mdr_tables, paste0(effect_var, "_MDR_final_fit.rds")))

  pred <- as_tibble(final_fit$model_output) %>%
    mutate(Date = pred_dates(., df_log), effect_var = effect_var) %>%
    relocate(Date, effect_var)
  coefs <- as_tibble(final_fit$smap_coefficients) %>%
    mutate(Date = coef_dates(., df_log), effect_var = effect_var) %>%
    relocate(Date, effect_var)

  save_csv(pred, file.path(paths$mdr_tables, paste0(effect_var, "_MDR_prediction.csv")))
  save_csv(coefs, file.path(paths$mdr_tables, paste0(effect_var, "_MDR_coefficients_raw.csv")))

  cmap <- coef_name_map(block_mvd, coefs)
  save_csv(cmap, file.path(paths$mdr_tables, paste0(effect_var, "_MDR_coefficient_column_map.csv")))
  if (any(!cmap$exists)) {
    warning(
      "Some expected coefficient columns were not found for ", effect_var, ": ",
      paste(cmap$coef_col[!cmap$exists], collapse = ", ")
    )
  }

  edge_coef <- purrr::map_dfr(seq_len(nrow(uic_for_block)), function(i) {
    link <- uic_for_block[i, ]
    adjusted_tp <- link$tp + config$smap_tp
    block_col <- paste0(link$cause_var, "_tp", adjusted_tp)
    coef_col <- cmap$coef_col[match(block_col, cmap$block_col)]

    if (length(coef_col) == 0 || is.na(coef_col) || !(coef_col %in% names(coefs))) {
      warning("Coefficient for ", link$edge, " not found. Expected block column: ", block_col)
      return(tibble())
    }

    tibble(
      Date = coefs$Date,
      effect_var = link$effect_var,
      cause_var = link$cause_var,
      edge = link$edge,
      lag_weeks = link$lag_weeks,
      uic_tp = link$tp,
      mdr_block_col = block_col,
      coef_col = coef_col,
      coef_value = coefs[[coef_col]]
    )
  })

  save_csv(edge_coef, file.path(paths$mdr_tables, paste0(effect_var, "_MDR_edge_coefficients_long.csv")))

  edge_summary <- edge_coef %>%
    group_by(.data$effect_var, .data$cause_var, .data$edge, .data$lag_weeks, .data$uic_tp, .data$mdr_block_col, .data$coef_col) %>%
    summarise(
      n_coef = sum(is.finite(.data$coef_value)),
      coef_mean = mean(.data$coef_value, na.rm = TRUE),
      coef_median = median(.data$coef_value, na.rm = TRUE),
      coef_sd = sd(.data$coef_value, na.rm = TRUE),
      coef_q25 = quantile(.data$coef_value, 0.25, na.rm = TRUE, names = FALSE),
      coef_q75 = quantile(.data$coef_value, 0.75, na.rm = TRUE, names = FALSE),
      prop_positive = mean(.data$coef_value > 0, na.rm = TRUE),
      prop_negative = mean(.data$coef_value < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      polarity = case_when(
        .data$coef_median > 0 ~ "positive",
        .data$coef_median < 0 ~ "negative",
        TRUE ~ "near_zero"
      ),
      effect_embedding_E = Ex,
      mdr_E = mvd_E,
      theta = best$theta[[1]],
      lambda = best$lambda[[1]],
      rho = best$rho[[1]],
      mae = best$mae[[1]],
      rmse = best$rmse[[1]],
      regularized = config$ridge_regularized_mdr
    )

  save_csv(edge_summary, file.path(paths$mdr_tables, paste0(effect_var, "_MDR_edge_summary.csv")))

  # Per-effect diagnostic plot: observed vs predicted.
  p_obs_pred <- ggplot(pred, aes(x = .data$obs, y = .data$pred)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.65, size = 1.7) +
    coord_equal() +
    labs(
      title = paste0("MDR S-map prediction: ", subtype_label(effect_var)),
      subtitle = paste0(
        "theta = ", signif(best$theta[[1]], 3),
        ", rho = ", signif(best$rho[[1]], 3),
        ", RMSE = ", signif(best$rmse[[1]], 3)
      ),
      x = "Observed log(x+1)",
      y = "Predicted log(x+1)"
    ) +
    theme_pub(13)
  ggsave(file.path(paths$mdr_fig, paste0(effect_var, "_MDR_observed_vs_predicted.tiff")), p_obs_pred,
    width = 7, height = 6, dpi = config$dpi, compression = "lzw"
  )

  # Per-effect coefficient time series.
  if (nrow(edge_coef) > 0) {
    p_coef_ts <- ggplot(edge_coef, aes(x = .data$Date, y = .data$coef_value)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_line(linewidth = 0.5, na.rm = TRUE) +
      facet_wrap(~edge, scales = "free_y", ncol = 1) +
      labs(
        title = paste0("Time-varying MDR S-map coefficients for ", subtype_label(effect_var)),
        x = NULL,
        y = "MDR S-map coefficient on log(x+1) scale"
      ) +
      theme_pub(13)
    ggsave(file.path(paths$mdr_fig, paste0(effect_var, "_MDR_edge_coefficients_timeseries.tiff")), p_coef_ts,
      width = 11, height = 4 + 2.5 * length(unique(edge_coef$edge)), dpi = config$dpi, compression = "lzw"
    )
  }

  list(
    effect_var = effect_var,
    Ex = Ex,
    block_mvd = block_mvd,
    dist_info = dist_info,
    param_res = param_res,
    best = best,
    pred = pred,
    coefs = coefs,
    coef_map = cmap,
    edge_coef = edge_coef,
    edge_summary = edge_summary
  )
}

mdr_results <- purrr::map(setNames(unique(selected_links$effect_var), unique(selected_links$effect_var)), function(effect_var) {
  run_mdr_for_effect(effect_var, selected_links, df_model, df_log, config, paths)
})
mdr_results <- purrr::compact(mdr_results)
save_rds(mdr_results, file.path(paths$mdr_tables, "MDR_results_all_effects.rds"))

edge_coef_all <- purrr::map_dfr(mdr_results, "edge_coef")
edge_summary_all <- purrr::map_dfr(mdr_results, "edge_summary")

# Join UIC information onto MDR summary.
edge_summary_all <- edge_summary_all %>%
  left_join(
    selected_links %>%
      select(
        effect_var, cause_var, E, tp, ete,
        q95_global, p_global, p_global_fdr
      ),
    by = c("effect_var", "cause_var")
  ) %>%
  rename(
    uic_E = E,
    uic_tp_selected = tp,
    uic_ete = ete,
    uic_q95_global = q95_global,
    uic_p_global = p_global,
    uic_p_global_fdr = p_global_fdr
  ) %>%
  arrange(.data$effect_var, .data$cause_var)

save_csv(edge_coef_all, file.path(paths$tables, "Figure5_MDR_coefficients_long.csv"))
save_csv(edge_summary_all, file.path(paths$tables, "Table_MDR_edge_summary.csv"))

# Compatibility export for downstream figure scripts. This does not overwrite old result folders.
coef_wide <- edge_coef_all %>%
  mutate(edge_col = paste0(.data$cause_var, "_cause_", .data$effect_var)) %>%
  group_by(.data$edge_col) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  select(row_id, edge_col, coef_value) %>%
  pivot_wider(names_from = edge_col, values_from = coef_value) %>%
  select(-row_id)
save_csv(coef_wide, file.path(paths$tables, "MDR_coefficients_all_wide_for_Figure5.csv"))

# -----------------------------
# 7. Manuscript-ready figures
# -----------------------------
fig5 <- edge_coef_all %>%
  mutate(edge = forcats::fct_reorder(.data$edge, .data$coef_value, .fun = median, na.rm = TRUE)) %>%
  ggplot(aes(x = .data$edge, y = .data$coef_value)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(alpha = 0.45, trim = FALSE, na.rm = TRUE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.75, na.rm = TRUE) +
  sinaplot::geom_sina(alpha = 0.35, size = 0.8, na.rm = TRUE) +
  coord_flip() +
  labs(
    title = "Figure 5. MDR S-map coefficients for UIC-selected subtype interactions",
    x = NULL,
    y = "MDR S-map coefficient on log(x+1) scale"
  ) +
  theme_pub(14)

ggsave(file.path(paths$fig, "Figure5_MDR_Smap_coefficients.tiff"), fig5,
  width = 10, height = max(5, 1.5 + length(unique(edge_coef_all$edge)) * 1.0),
  dpi = config$dpi, compression = "lzw"
)

# Figure 4 network summary.
coords <- tibble(
  node = c("B", "A_H1N1", "A_H3N2"),
  x = c(0, -1, 1),
  y = c(1, 0, 0)
)

edges_net <- edge_summary_all %>%
  left_join(coords, by = c("cause_var" = "node")) %>%
  rename(x0 = x, y0 = y) %>%
  left_join(coords, by = c("effect_var" = "node")) %>%
  rename(x1 = x, y1 = y) %>%
  mutate(
    label = paste0("lag ", .data$lag_weeks, " wk\n", sprintf("median %.3f", .data$coef_median)),
    abs_coef = abs(.data$coef_median)
  )
abs_rng <- range(edges_net$abs_coef, na.rm = TRUE)
if (!all(is.finite(abs_rng)) || abs_rng[1] == abs_rng[2]) {
  edges_net$line_width <- 1.5
} else {
  edges_net$line_width <- scales::rescale(edges_net$abs_coef, to = c(0.5, 2.6), from = abs_rng)
}

save_csv(edges_net, file.path(paths$tables, "Figure4_network_edges_MDR.csv"))

fig4 <- ggplot() +
  geom_segment(
    data = edges_net,
    aes(
      x = .data$x0, y = .data$y0, xend = .data$x1, yend = .data$y1,
      linewidth = .data$line_width, linetype = .data$polarity
    ),
    arrow = grid::arrow(length = unit(0.22, "cm"), type = "closed"),
    alpha = 0.8
  ) +
  geom_point(data = coords, aes(x = .data$x, y = .data$y), size = 16, shape = 21, fill = "white") +
  geom_text(data = coords, aes(x = .data$x, y = .data$y, label = subtype_label(.data$node)), size = 5, fontface = "bold") +
  geom_text(data = edges_net, aes(x = (.data$x0 + .data$x1) / 2, y = (.data$y0 + .data$y1) / 2 + 0.08, label = .data$label), size = 3.4) +
  scale_linewidth_identity() +
  coord_equal(xlim = c(-1.6, 1.6), ylim = c(-0.25, 1.25)) +
  labs(
    title = "Figure 4. Directional subtype interaction network from UIC and MDR S-map",
    linetype = "median coefficient"
  ) +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")

ggsave(file.path(paths$fig, "Figure4_UIC_MDR_network.tiff"), fig4,
  width = 8, height = 6, dpi = config$dpi, compression = "lzw"
)

# -----------------------------
# 8. Manuscript-ready tables and text snippets
# -----------------------------
table_main <- edge_summary_all %>%
  transmute(
    Cause = subtype_label(.data$cause_var),
    Effect = subtype_label(.data$effect_var),
    `UIC lag, weeks` = .data$lag_weeks,
    `UIC E` = .data$uic_E,
    `UIC ETE` = round(.data$uic_ete, 4),
    `UIC FDR p` = signif(.data$uic_p_global_fdr, 3),
    `MDR theta` = signif(.data$theta, 3),
    `MDR rho` = signif(.data$rho, 3),
    `MDR RMSE` = signif(.data$rmse, 3),
    `Mean MDR coef` = round(.data$coef_mean, 4),
    `Median MDR coef` = round(.data$coef_median, 4),
    `IQR` = paste0(round(.data$coef_q25, 4), " to ", round(.data$coef_q75, 4)),
    `Positive proportion` = round(.data$prop_positive, 3),
    Polarity = .data$polarity
  )

save_csv(table_main, file.path(paths$tables, "Table1_UIC_MDR_main_results.csv"))

methods_text <- c(
  "Suggested Methods replacement for the S-map paragraph:",
  "",
  paste0(
    "Interaction strength and polarity were estimated using multiview-distance S-map (MDR S-map). ",
    "For each effect subtype, UIC was first performed across the remaining subtypes over lags of ",
    paste(range(config$tp_range), collapse = " to "),
    ". UIC links were retained for MDR S-map only when the observed effective transfer entropy exceeded ",
    "the 95th percentile of the 52-week seasonal surrogate max-statistic distribution and the empirical ",
    "surrogate p-value remained significant after Benjamini-Hochberg false-discovery-rate correction. ",
    "The retained UIC links were then used to build the MDR block with make_block_mvd(). Because the MDR S-map ",
    "was fitted as a one-step map, UIC tp values were adjusted by +1 in the block construction so that a selected ",
    "UIC lag of k weeks was represented as the local coefficient of the cause at t-k+1 for the prediction of the ",
    "effect at t+1. Multiview distances were calculated using compute_mvd() with E = ", config$mdr_E,
    ", n_ssr = ", config$n_ssr, ", and k = ", config$k,
    ". The S-map weighting parameter theta was selected by minimizing RMSE over the prespecified grid. ",
    "MDR S-map coefficients were interpreted on the log(x+1)-transformed incidence scale."
  ),
  "",
  "Suggested Figure 5 legend:",
  "Figure 5. Time-varying MDR S-map coefficients for UIC-selected subtype interactions. Coefficients quantify the local effect of the cause subtype on the subsequent dynamics of the effect subtype at the UIC-selected lag, on the log(x+1)-transformed scale. Violin plots show the full empirical distribution, boxes show the interquartile range and median, and the dashed horizontal line indicates zero.",
  "",
  "Suggested Results table is exported as tables/Table1_UIC_MDR_main_results.csv. Replace numerical values in the manuscript after running the script on the final dataset."
)
writeLines(methods_text, file.path(paths$manuscript, "MDR_methods_results_figure_legend_text.txt"))

# -----------------------------
# 9. Session information
# -----------------------------
sink(file.path(paths$root, "sessionInfo.txt"))
cat("MDR package used: ", mdr_pkg, "\n")
cat("MDR package version: ", as.character(utils::packageVersion(mdr_pkg)), "\n")
cat("rUIC version: ", as.character(utils::packageVersion("rUIC")), "\n")
cat("rEDM version: ", as.character(utils::packageVersion("rEDM")), "\n")
cat("\nConfiguration:\n")
print(config)
cat("\nSession info:\n")
print(sessionInfo())
sink()

message("Done. Results written to: ", normalizePath(paths$root, winslash = "/", mustWork = FALSE))
