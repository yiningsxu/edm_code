# flu_subtype_mdr_smap_uic_optimal_noFDR.R
# Primary EDM pipeline for influenza subtype interaction analysis using UIC + MDR S-map.
# Variant: NO Benjamini-Hochberg/FDR correction is applied to UIC p-values.
#
# Purpose:
#   1) Detect directional, lagged subtype interactions by UIC with seasonal surrogate correction.
#   2) Use max-statistic seasonal-surrogate-significant UIC links to construct multiview-distance blocks.
#   3) Estimate interaction polarity/strength with MDR S-map as the primary S-map analysis.
#   4) Export manuscript-ready tables, Figure 3/4/5 inputs, and figure files.
#
# Expected input:
#   data/FluSub_jp/FluSub_11to19_jp_per_20240925.csv
#   Required columns: year, week, B, A_H1N1, A_H3N2
#
# Notes:
#   - Run with Rscript or source() from anywhere inside the edm_code workspace.
#   - For a quick dry run, use environment variables such as
#     FLU_SUBTYPE_NUM_SURR=10 FLU_SUBTYPE_UIC_INTERNAL_NUM_SURR=10 FLU_SUBTYPE_N_SSR=50 FLU_SUBTYPE_E_RANGE=1:5.
#   - MDR distance and S-map functions are taken from either macam or macamts.
#   - This variant uses rUIC::uic.optimal() for UIC, not rUIC::uic() with lag-wise max-ETE E selection.
#   - This no-FDR variant retains UIC links when p_global < alpha and ETE > q95_global.

rm(list = ls())

# -----------------------------
# 0. Bootstrap and packages
# -----------------------------
current_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- args[startsWith(args, "--file=")][1]
    if (!is.na(file_arg)) {
        file_arg <- sub("^--file=", "", file_arg)
        if (file.exists(file_arg)) {
            return(dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE)))
        }
    }

    for (frame in rev(sys.frames())) {
        ofile <- frame$ofile
        if (!is.null(ofile) && nzchar(ofile) && file.exists(ofile)) {
            return(dirname(normalizePath(ofile, winslash = "/", mustWork = TRUE)))
        }
    }

    NA_character_
}

# ---- edm_code bootstrap ----
source_edm_bootstrap <- function(start_dirs = c(current_script_dir(), getwd())) {
    start_dirs <- unique(start_dirs[!is.na(start_dirs) & nzchar(start_dirs)])

    for (start_dir in start_dirs) {
        current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
        repeat {
            bootstrap <- file.path(current, "R", "bootstrap.R")
            if (file.exists(bootstrap)) {
                source(bootstrap)
                return(source_edm_paths(current))
            }
            parent <- dirname(current)
            if (identical(parent, current)) {
                break
            }
            current <- parent
        }
    }

    stop("Could not find edm_code/R/bootstrap.R. Run this script from inside edm_code.", call. = FALSE)
}
source_edm_bootstrap()
setwd(workspace_root())
rm(source_edm_bootstrap, current_script_dir)

required_packages <- c(
    "dplyr", # データ加工
    "tidyr", # データの縦横変換
    "purrr", # 関数型処理、反復処理
    "tibble", # データフレームの拡張
    "readr", # csvなどの読み込み
    "stringr", # 文字列処理
    "forcats", # カテゴリ変数の処理
    "scales", # ggplot の軸・ラベル調整
    "lubridate", # 日付・時刻の操作
    "ISOweek", # ISO週番号の計算
    "ggplot2", # グラフ作成
    "cowplot", # 複数図の結合
    "ggforce", # sina plot layer
    "rUIC", # Unified Information Criterion 関連
    "rEDM", # Empirical Dynamic Modeling
    "grid" # 図の細かい制御
)
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0) {
    stop(
        "Missing required R package(s): ",
        paste(missing_packages, collapse = ", "),
        ". Install them before running this pipeline.",
        call. = FALSE
    )
}
invisible(lapply(required_packages, library, character.only = TRUE))
# Keep dplyr at the front of the search path so unqualified verbs retain
# dplyr semantics even if an attached dependency masks them.
if ("package:dplyr" %in% search()) {
    detach("package:dplyr", unload = FALSE, character.only = TRUE)
}
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
rm(required_packages, missing_packages)

message("rEDM version: ", as.character(utils::packageVersion("rEDM")))
message("rUIC version: ", as.character(utils::packageVersion("rUIC")))

resolve_mdr_pkg <- function() {
    candidates <- c("macam", "macamts")
    for (pkg in candidates) {
        has_namespace <- requireNamespace(pkg, quietly = TRUE)
        has_functions <- has_namespace &&
            all(c("compute_mvd", "s_map_mdr") %in% getNamespaceExports(pkg))
        if (has_functions) {
            return(pkg)
        }
    }
    stop(
        "Neither 'macam' nor 'macamts' with compute_mvd() and s_map_mdr() is installed. ",
        "Install one of them, e.g. remotes::install_github('ong8181/macam').",
        call. = FALSE
    )
}

mdr_pkg <- resolve_mdr_pkg()
message("Using MDR package: ", mdr_pkg, " ", as.character(utils::packageVersion(mdr_pkg)))

mdr_fun <- function(fun) {
    getExportedValue(mdr_pkg, fun)
}

# Local make_block_mvd wrapper.
# The installed macam/macamts versions in this workspace do not expose tp_adjust,
# but tp_adjust is needed to align a one-step MDR S-map target with UIC lags.
make_block_mvd_compat <- function(block,
                                  uic_res,
                                  effect_var,
                                  max_lag,
                                  cause_var_colname = "cause_var",
                                  include_var = "strongest_only",
                                  p_threshold = 0.05,
                                  tp_adjust = 0,
                                  sort_tp = TRUE,
                                  silent = FALSE,
                                  ...) {
    dots <- list(...)
    if (length(dots) > 0) {
        warning("Unused make_block_mvd_compat() argument(s): ", paste(names(dots), collapse = ", "))
    }

    block <- as.data.frame(block)
    uic_res <- as.data.frame(uic_res)
    x_names <- colnames(block)

    if (is.numeric(effect_var)) {
        effect_var <- x_names[effect_var]
    }
    if (!(effect_var %in% x_names)) {
        stop("No effect_var in the block: ", effect_var, call. = FALSE)
    }
    if (!all(unique(x_names) == x_names)) {
        stop("block must have unique column names.", call. = FALSE)
    }
    if (!include_var %in% c("all_significant", "strongest_only", "tp0_only")) {
        stop("include_var must be all_significant, strongest_only, or tp0_only.", call. = FALSE)
    }
    if (!is.numeric(max_lag) || length(max_lag) != 1 || max_lag < 1) {
        stop("max_lag must be a single number >= 1.", call. = FALSE)
    }

    required_uic_cols <- c(cause_var_colname, "tp", "te", "pval")
    missing_uic_cols <- setdiff(required_uic_cols, colnames(uic_res))
    if (length(missing_uic_cols) > 0) {
        stop("uic_res is missing required column(s): ", paste(missing_uic_cols, collapse = ", "), call. = FALSE)
    }

    max_lag <- as.integer(max_lag)
    tp_adjust <- as.integer(tp_adjust)

    if (include_var == "tp0_only") {
        block_mvd <- data.frame(block[[effect_var]])
        colnames(block_mvd) <- sprintf("%s_tp0", effect_var)
    } else {
        block_mvd <- data.frame(rEDM::make_block(block[[effect_var]], max_lag = max_lag)[, -1, drop = FALSE])
        colnames(block_mvd) <- sprintf("%s_tp%s", effect_var, 0:(-(max_lag - 1)))
    }

    if (!silent) {
        message(sprintf(
            "UIC rows with tp <= 0 and pval <= %s are kept; cause tp is adjusted by %+d for MDR.",
            p_threshold,
            tp_adjust
        ))
    }

    uic_res <- uic_res[is.finite(uic_res$pval) & uic_res$pval <= p_threshold & uic_res$tp <= 0, , drop = FALSE]
    if (nrow(uic_res) < 1) {
        stop("No significant causal variables were detected. Please use the univariate S-map.", call. = FALSE)
    }

    if (sort_tp) {
        uic_res <- uic_res %>%
            dplyr::group_by(.data[[cause_var_colname]]) %>%
            dplyr::arrange(dplyr::desc(.data$tp), .by_group = TRUE) %>%
            dplyr::ungroup() %>%
            as.data.frame()
    }

    if (include_var == "strongest_only") {
        uic_res <- uic_res %>%
            dplyr::group_by(.data[[cause_var_colname]]) %>%
            dplyr::arrange(dplyr::desc(.data$te), .by_group = TRUE) %>%
            dplyr::slice(1) %>%
            dplyr::ungroup() %>%
            as.data.frame()
    }

    if (include_var == "tp0_only") {
        for (cause_i in unique(uic_res[[cause_var_colname]])) {
            if (!(cause_i %in% x_names)) {
                stop("UIC cause variable is not in block: ", cause_i, call. = FALSE)
            }
            block_mvd[[sprintf("%s_tp0", cause_i)]] <- block[[cause_i]]
        }
        return(block_mvd)
    }

    for (i in seq_len(nrow(uic_res))) {
        cause_i <- uic_res[[cause_var_colname]][i]
        if (!(cause_i %in% x_names)) {
            stop("UIC cause variable is not in block: ", cause_i, call. = FALSE)
        }

        adjusted_tp <- as.integer(uic_res$tp[i]) + tp_adjust
        if (adjusted_tp > 0) {
            stop(
                "Adjusted cause tp became positive for ", cause_i, " -> ", effect_var,
                " (UIC tp = ", uic_res$tp[i], ", tp_adjust = ", tp_adjust, ").",
                call. = FALSE
            )
        }

        block_new <- if (adjusted_tp == 0) {
            block[[cause_i]]
        } else {
            dplyr::lag(block[[cause_i]], n = abs(adjusted_tp))
        }
        block_mvd[[sprintf("%s_tp%s", cause_i, adjusted_tp)]] <- block_new
    }

    block_mvd
}

call_mdr_function <- function(fun, ...) {
    f <- mdr_fun(fun) # MDR パッケージから関数を取得する
    args <- list(...) # 渡された引数をリスト化する
    formal_names <- names(formals(f)) # 関数の引数名を取得する

    dropped <- setdiff(names(args), formal_names)
    if (length(dropped) > 0) {
        warning(
            fun, "() does not support argument(s): ",
            paste(dropped, collapse = ", "),
            ". These were dropped."
        )
    }


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

analysis_start_time <- Sys.time()

log_msg <- function(..., level = "INFO") {
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    print(paste0("[", stamp, "] [", level, "] ", paste0(..., collapse = "")))
    flush.console()
}

log_section <- function(title) {
    print("")
    log_msg("==== ", title, " ====")
}

format_elapsed <- function(start_time) {
    secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (!is.finite(secs)) {
        return("unknown")
    }
    if (secs < 60) {
        return(paste0(round(secs, 1), " sec"))
    }
    if (secs < 3600) {
        return(paste0(round(secs / 60, 1), " min"))
    }
    paste0(round(secs / 3600, 2), " hr")
}

progress_points <- function(total, n = 10) {
    if (total <= 0) {
        return(integer(0))
    }
    sort(unique(pmax(1L, pmin(total, round(seq(1, total, length.out = min(n, total)))))))
}

# -----------------------------
# 1. Configuration
# -----------------------------
log_section("1. Configuration")

env_chr <- function(name, default) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (is.na(value) || !nzchar(value)) {
        return(default)
    }
    value
}

env_int <- function(name, default) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (is.na(value) || !nzchar(value)) {
        return(default)
    }
    out <- suppressWarnings(as.integer(value))
    if (is.na(out)) {
        stop("Environment variable ", name, " must be an integer.", call. = FALSE)
    }
    out
}

env_num <- function(name, default) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (is.na(value) || !nzchar(value)) {
        return(default)
    }
    out <- suppressWarnings(as.numeric(value))
    if (is.na(out)) {
        stop("Environment variable ", name, " must be numeric.", call. = FALSE)
    }
    out
}

env_int_nullable <- function(name, default = NULL) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (is.na(value) || !nzchar(value)) {
        return(default)
    }
    out <- suppressWarnings(as.integer(value))
    if (is.na(out)) {
        stop("Environment variable ", name, " must be an integer.", call. = FALSE)
    }
    out
}

env_int_vector <- function(name, default) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (is.na(value) || !nzchar(value)) {
        return(default)
    }
    value <- gsub("\\s+", "", value)
    if (grepl("^-?[0-9]+:-?[0-9]+$", value)) {
        bounds <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
        return(seq(bounds[1], bounds[2]))
    }
    out <- suppressWarnings(as.integer(strsplit(value, ",", fixed = TRUE)[[1]]))
    if (anyNA(out)) {
        stop("Environment variable ", name, " must be an integer range such as 1:5 or a comma-separated integer list.", call. = FALSE)
    }
    out
}

env_num_vector <- function(name, default) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (is.na(value) || !nzchar(value)) {
        return(default)
    }
    out <- suppressWarnings(as.numeric(strsplit(gsub("\\s+", "", value), ",", fixed = TRUE)[[1]]))
    if (anyNA(out)) {
        stop("Environment variable ", name, " must be a comma-separated numeric list.", call. = FALSE)
    }
    out
}

env_logical <- function(name, default) {
    value <- tolower(Sys.getenv(name, unset = NA_character_))
    if (is.na(value) || !nzchar(value)) {
        return(default)
    }
    if (value %in% c("true", "t", "1", "yes", "y")) {
        return(TRUE)
    }
    if (value %in% c("false", "f", "0", "no", "n")) {
        return(FALSE)
    }
    stop("Environment variable ", name, " must be TRUE or FALSE.", call. = FALSE)
}

config <- list(
    data_file = env_chr("FLU_SUBTYPE_DATA_FILE", "data/FluSub_jp/FluSub_11to19_jp_per_20240925.csv"),
    out_dir = env_chr(
        "FLU_SUBTYPE_OUT_DIR",
        file.path("result", "FluSub_JP", paste0(format(Sys.Date(), "%Y%m%d"), "_UICoptimal_MDR_noFDR"))
    ),
    subtype_vars = c("B", "A_H1N1", "A_H3N2"),

    # UIC settings. Primary analysis excludes tp = 0 to avoid contemporaneous seasonal synchrony.
    E_range = env_int_vector("FLU_SUBTYPE_E_RANGE", 1:20), # 埋め込み次元またはラグ数候補の範囲
    tp_range = env_int_vector("FLU_SUBTYPE_TP_RANGE", -12:-1), # 予測ラグの範囲
    tau = env_int("FLU_SUBTYPE_TAU", 1), # 予測ステップ数（デフォルトは1, 1週間刻みでラグを取る）
    alpha = env_num("FLU_SUBTYPE_ALPHA", 0.05), # 有意水準
    num_surr = env_int("FLU_SUBTYPE_NUM_SURR", 2000), # 季節サロゲートデータの数
    uic_internal_num_surr = env_int("FLU_SUBTYPE_UIC_INTERNAL_NUM_SURR", 1000), # rUIC::uic.optimal() 内部p値用。dry runでは小さくしてよい
    season_period = env_int("FLU_SUBTYPE_SEASON_PERIOD", 52), # 季節周期
    random_seed = env_int("FLU_SUBTYPE_RANDOM_SEED", 1234), # 乱数シード
    save_intermediate = env_logical("FLU_SUBTYPE_SAVE_INTERMEDIATE", TRUE), # 解析途中のcheckpoint保存
    save_each_surrogate = env_logical("FLU_SUBTYPE_SAVE_EACH_SURROGATE", TRUE), # サロゲートごとのUIC結果を個別保存
    checkpoint_every = env_int("FLU_SUBTYPE_CHECKPOINT_EVERY", 25), # サロゲート行列などの途中保存間隔

    # MDR S-map settings.
    smap_tp = env_int("FLU_SUBTYPE_SMAP_TP", 1), # S-map の予測ターゲットとなる時間ステップ
    mdr_include_var = "strongest_only", # 各効果変数に対して、UICで最も強い関係を示した原因変数だけを MDR S-map に入れる設定
    mdr_E = env_int("FLU_SUBTYPE_MDR_E", 3), # MDR S-map に使用する埋め込み次元
    n_ssr = env_int("FLU_SUBTYPE_N_SSR", 2000), # random multiview embeddingの候補数
    k = env_int_nullable("FLU_SUBTYPE_K", NULL), # if NULL, floor(sqrt(n_ssr)) is used.
    theta_grid = env_num_vector(
        "FLU_SUBTYPE_THETA_GRID",
        c(0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
    ), # S-mapでは、theta が大きいほど近い状態にある点を重視する

    # 正則化付き MDR S-map を使うかどうか。感度分析として、TRUE のままで解析を進めることも可能
    ridge_regularized_mdr = env_logical("FLU_SUBTYPE_RIDGE_REGULARIZED_MDR", FALSE), # 主解析は使わない
    lambda_grid = env_num_vector("FLU_SUBTYPE_LAMBDA_GRID", c(0)), # 正則化の強さのグリッド値。0は正則化なし
    alpha_glmnet = 0, # alpha = 0：Ridge, alpha = 1：Lasso, 0 < alpha < 1：Elastic Net

    # Figure export.
    dpi = env_int("FLU_SUBTYPE_DPI", 300),
    fig_width = 12,
    fig_height = 8
)
if (is.null(config$k)) config$k <- floor(sqrt(config$n_ssr))
if (config$num_surr < 1) stop("config$num_surr must be >= 1.", call. = FALSE)
if (config$uic_internal_num_surr < 1) stop("config$uic_internal_num_surr must be >= 1.", call. = FALSE)
if (config$n_ssr < 1) stop("config$n_ssr must be >= 1.", call. = FALSE)
if (config$k < 1) stop("config$k must be >= 1.", call. = FALSE)
if (config$alpha <= 0 || config$alpha > 1) stop("config$alpha must be in (0, 1].", call. = FALSE)
if (config$checkpoint_every < 1) stop("config$checkpoint_every must be >= 1.", call. = FALSE)
log_msg("Input file: ", config$data_file)
log_msg("Subtypes: ", paste(config$subtype_vars, collapse = ", "))
log_msg(
    "UIC settings: E = ", paste(range(config$E_range), collapse = "-"),
    ", tp = ", paste(range(config$tp_range), collapse = " to "),
    ", seasonal surrogates = ", config$num_surr,
    ", rUIC uic.optimal internal surrogates = ", config$uic_internal_num_surr
)
log_msg(
    "Intermediate saving: save_intermediate = ", config$save_intermediate,
    ", save_each_surrogate = ", config$save_each_surrogate,
    ", checkpoint_every = ", config$checkpoint_every
)
log_msg(
    "MDR settings: mdr_E = ", config$mdr_E,
    ", n_ssr = ", config$n_ssr,
    ", k = ", config$k,
    ", theta candidates = ", length(config$theta_grid)
)

# -----------------------------
# 2. Output paths
# -----------------------------
log_section("2. Output paths")
safe_dir <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
}

paths <- list(
    root = safe_dir(config$out_dir),
    tables = safe_dir(file.path(config$out_dir, "tables")),
    uic_tables = safe_dir(file.path(config$out_dir, "tables", "uic")),
    uic_surrogate_tables = safe_dir(file.path(config$out_dir, "tables", "uic", "surrogate_runs")),
    mdr_tables = safe_dir(file.path(config$out_dir, "tables", "mdr_smap")),
    checkpoints = safe_dir(file.path(config$out_dir, "checkpoints")),
    uic_checkpoints = safe_dir(file.path(config$out_dir, "checkpoints", "uic")),
    mdr_checkpoints = safe_dir(file.path(config$out_dir, "checkpoints", "mdr_smap")),
    progress = safe_dir(file.path(config$out_dir, "progress")),
    fig = safe_dir(file.path(config$out_dir, "figures")),
    uic_fig = safe_dir(file.path(config$out_dir, "figures", "uic")),
    mdr_fig = safe_dir(file.path(config$out_dir, "figures", "mdr_smap")),
    manuscript = safe_dir(file.path(config$out_dir, "manuscript_ready"))
)
log_msg("Output root: ", normalizePath(paths$root, winslash = "/", mustWork = FALSE))
log_msg("Tables: ", normalizePath(paths$tables, winslash = "/", mustWork = FALSE))
log_msg("Figures: ", normalizePath(paths$fig, winslash = "/", mustWork = FALSE))

save_csv <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(x, path, na = "")
    invisible(path)
}

save_rds <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(x, path)
    invisible(path)
}

append_csv <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    append <- file.exists(path)
    readr::write_csv(x, path, na = "", append = append, col_names = !append)
    invisible(path)
}

reset_output_file <- function(path) {
    if (file.exists(path)) {
        file.remove(path)
    }
    invisible(path)
}

config_tbl <- function(config) {
    tibble(
        name = names(config),
        value = vapply(config, function(x) paste(as.character(x), collapse = ","), character(1))
    )
}

log_progress <- function(stage, tag = NA_character_, detail = NA_character_, path = NA_character_) {
    append_csv(
        tibble(
            timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            elapsed = format_elapsed(analysis_start_time),
            stage = stage,
            tag = tag,
            detail = detail,
            path = path
        ),
        file.path(paths$progress, "progress_log.csv")
    )
}

save_status <- function(stage, status, detail = NA_character_, path = NA_character_) {
    status_tbl <- tibble(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        elapsed = format_elapsed(analysis_start_time),
        stage = stage,
        status = status,
        detail = detail,
        path = path
    )
    save_csv(status_tbl, file.path(paths$progress, "analysis_status_latest.csv"))
    log_progress(stage, tag = status, detail = detail, path = path)
    invisible(status_tbl)
}

write_session_info <- function(status = "completed", notes = character()) {
    lines <- capture.output({
        cat("Analysis status: ", status, "\n", sep = "")
        if (length(notes) > 0) {
            cat("\nNotes:\n")
            cat(paste0("- ", notes, collapse = "\n"), "\n")
        }
        cat("\nMDR package used: ", mdr_pkg, "\n", sep = "")
        cat("MDR package version: ", as.character(utils::packageVersion(mdr_pkg)), "\n", sep = "")
        cat("rUIC version: ", as.character(utils::packageVersion("rUIC")), "\n", sep = "")
        cat("rEDM version: ", as.character(utils::packageVersion("rEDM")), "\n", sep = "")
        cat("\nConfiguration:\n")
        print(config)
        cat("\nSession info:\n")
        print(sessionInfo())
    })
    writeLines(lines, file.path(paths$root, "sessionInfo.txt"))
    invisible(file.path(paths$root, "sessionInfo.txt"))
}

reset_output_file(file.path(paths$progress, "progress_log.csv"))
save_csv(config_tbl(config), file.path(paths$tables, "run_config.csv"))
save_rds(list(config = config, paths = paths), file.path(paths$root, "run_config_and_paths.rds"))
save_status(
    stage = "configuration",
    status = "saved",
    detail = "Run configuration and output paths were saved.",
    path = file.path(paths$root, "run_config_and_paths.rds")
)

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
log_section("3. Data preparation")
read_prepare_flu <- function(data_file, subtype_vars) {
    log_msg("Reading influenza subtype data: ", data_file)
    if (!file.exists(data_file)) {
        stop("Data file not found: ", data_file, call. = FALSE)
    }

    raw <- read.csv(data_file, check.names = FALSE)
    log_msg("Raw data loaded: ", nrow(raw), " rows x ", ncol(raw), " columns")

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
    log_msg(
        "Date range after ISO-week conversion: ",
        min(raw$Date, na.rm = TRUE), " to ", max(raw$Date, na.rm = TRUE)
    )

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
        if (any(x < 0)) {
            stop("Subtype variable contains negative values: ", v, call. = FALSE)
        }

        out[[v]] <- log1p(x)
    }

    keep <- c(
        TRUE,
        vapply(out[subtype_vars], function(x) sum(x, na.rm = TRUE) > 0, logical(1))
    )

    out <- out[, keep, drop = FALSE]
    dropped_vars <- setdiff(subtype_vars, names(out))
    if (length(dropped_vars) > 0) {
        log_msg("Dropped all-zero subtype columns: ", paste(dropped_vars, collapse = ", "), level = "WARN")
    }

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
log_msg(
    "Prepared model data: ", nrow(df_model), " time points; variables = ",
    paste(vars, collapse = ", ")
)
log_msg("Saved prepared time series: ", file.path(paths$tables, "prepared_log1p_timeseries.csv"))
save_status(
    stage = "data_preparation",
    status = "completed",
    detail = paste0(nrow(df_model), " time points; variables: ", paste(vars, collapse = ", ")),
    path = file.path(paths$tables, "prepared_log1p_timeseries.csv")
)

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
log_section("5. UIC with seasonal surrogate correction")

run_uic_optimal <- function(block, lib_var, tar_var, E, tau, tp, alpha, num_surr = 1000) {
    # Use rUIC::uic.optimal() so that the embedding dimension is selected by
    # the rUIC wrapper before lag-wise UIC values are returned. This avoids the
    # previous behavior of computing all E x tp combinations with rUIC::uic()
    # and then selecting the E with the maximum ETE separately for each tp.
    f <- rUIC::uic.optimal
    formal_names <- names(formals(f))

    args <- list(as.data.frame(block))
    named_args <- list(
        lib_var = lib_var,
        tar_var = tar_var,
        E = E,
        tau = tau,
        tp = tp,
        alpha = alpha,
        num_surr = num_surr,
        sequential_test = FALSE
    )
    named_args <- named_args[names(named_args) %in% formal_names]

    out <- do.call(f, c(args, named_args)) %>%
        as_tibble() %>%
        dplyr::filter(.data$tp %in% .env$tp) %>%
        dplyr::arrange(.data$tp)

    if (!("tp" %in% names(out))) {
        stop("rUIC::uic.optimal() output does not include 'tp'.", call. = FALSE)
    }
    if (!("ete" %in% names(out))) {
        stop("rUIC::uic.optimal() output does not include 'ete'.", call. = FALSE)
    }
    if (!("te" %in% names(out))) {
        out <- out %>% dplyr::mutate(te = .data$ete)
    }
    if (!("pval" %in% names(out))) {
        out <- out %>% dplyr::mutate(pval = NA_real_)
    }
    if (!("E" %in% names(out))) {
        out <- out %>% dplyr::mutate(E = NA_integer_)
    }

    expected_tp <- sort(unique(tp))
    observed_tp <- sort(unique(out$tp))
    missing_tp <- setdiff(expected_tp, observed_tp)
    if (length(missing_tp) > 0) {
        stop(
            "rUIC::uic.optimal() did not return all requested tp values. Missing: ",
            paste(missing_tp, collapse = ", "),
            call. = FALSE
        )
    }

    duplicate_tp <- unique(out$tp[duplicated(out$tp)])
    if (length(duplicate_tp) > 0) {
        stop(
            "rUIC::uic.optimal() returned multiple rows for tp value(s): ",
            paste(duplicate_tp, collapse = ", "),
            call. = FALSE
        )
    }

    out
}

run_uic_pair <- function(df_model, effect_var, cause_var, config, paths, pair_id = NULL, total_pairs = NULL) {
    pair_start <- Sys.time()
    pair_prefix <- if (!is.null(pair_id) && !is.null(total_pairs)) {
        paste0("pair ", pair_id, "/", total_pairs, ": ")
    } else {
        ""
    }
    log_msg("UIC started: ", pair_prefix, cause_var, " -> ", effect_var)
    tag <- paste0(cause_var, "_to_", effect_var)
    pair_checkpoint_dir <- safe_dir(file.path(paths$uic_checkpoints, tag))
    surrogate_run_dir <- safe_dir(file.path(paths$uic_surrogate_tables, tag))
    surrogate_file_id <- function(i) {
        sprintf(paste0("%0", nchar(as.character(config$num_surr)), "d"), i)
    }
    log_progress(
        stage = "uic_pair_started",
        tag = tag,
        detail = paste0(pair_prefix, cause_var, " -> ", effect_var)
    )

    obs <- run_uic_optimal(
        block = df_model,
        lib_var = effect_var,
        tar_var = cause_var,
        E = config$E_range,
        tau = config$tau,
        tp = config$tp_range,
        alpha = config$alpha,
        num_surr = config$uic_internal_num_surr
    )
    log_msg(
        "Observed UIC completed for ", cause_var, " -> ", effect_var,
        ": ", nrow(obs), " tp rows retained"
    )

    if (!("ete" %in% names(obs))) {
        stop("rUIC::uic.optimal() output does not include 'ete'.", call. = FALSE)
    }
    if (!("te" %in% names(obs))) {
        obs <- obs %>% mutate(te = .data$ete)
    }
    if (!("pval" %in% names(obs))) {
        obs <- obs %>% mutate(pval = NA_real_)
    }
    if (isTRUE(config$save_intermediate)) {
        observed_path <- file.path(pair_checkpoint_dir, paste0(tag, "_observed_uic.csv"))
        save_csv(obs, observed_path)
        log_msg("Saved observed UIC checkpoint for ", tag)
        log_progress("uic_observed_saved", tag, paste0(nrow(obs), " tp rows"), observed_path)
    }
    obs_best <- obs %>%
        arrange(desc(.data$ete)) %>%
        slice(1)
    selected_E_values <- unique(stats::na.omit(obs$E))
    selected_E_msg <- if (length(selected_E_values) == 1) {
        as.character(selected_E_values[[1]])
    } else {
        paste(selected_E_values, collapse = ", ")
    }
    log_msg(
        "uic.optimal selected E = ", selected_E_msg,
        "; observed best lag: tp = ", obs_best$tp[[1]],
        ", ete = ", signif(obs_best$ete[[1]], 4)
    )

    log_msg(
        "Creating ", config$num_surr, " seasonal surrogates for effect variable ",
        effect_var, " (period = ", config$season_period, ")"
    )
    effect_surr <- make_seasonal_surrogates(
        df_model[[effect_var]],
        num_surr = config$num_surr,
        period = config$season_period,
        seed = config$random_seed
    )
    log_msg("Seasonal surrogate matrix created: ", nrow(effect_surr), " x ", ncol(effect_surr))
    if (isTRUE(config$save_intermediate)) {
        effect_surr_tbl <- tibble(time = seq_len(nrow(effect_surr))) %>%
            bind_cols(as_tibble(effect_surr, .name_repair = "minimal"))
        surrogate_series_path <- file.path(pair_checkpoint_dir, paste0(tag, "_seasonal_surrogate_series.csv"))
        save_csv(effect_surr_tbl, surrogate_series_path)
        log_msg("Saved seasonal surrogate series checkpoint for ", tag)
        log_progress("uic_surrogate_series_saved", tag, paste0(config$num_surr, " surrogate series"), surrogate_series_path)
    }

    surr_ete <- matrix(NA_real_, nrow = nrow(obs), ncol = config$num_surr)
    rownames(surr_ete) <- paste0("tp_", obs$tp)
    colnames(surr_ete) <- paste0("surr_", seq_len(config$num_surr))
    surrogate_errors <- list()
    surrogate_success <- 0L

    log_msg("Running surrogate uic.optimal loop for ", cause_var, " -> ", effect_var)
    surr_progress <- progress_points(config$num_surr, n = 10)
    for (i in seq_len(config$num_surr)) {
        tmp <- data.frame(effect = effect_surr[, i], cause = df_model[[cause_var]])
        sres <- tryCatch(
            {
                run_uic_optimal(
                    block = tmp,
                    lib_var = "effect",
                    tar_var = "cause",
                    E = config$E_range,
                    tau = config$tau,
                    tp = config$tp_range,
                    alpha = config$alpha,
                    num_surr = config$uic_internal_num_surr
                ) %>%
                    dplyr::select(any_of(c("tp", "E", "te", "ete", "pval")))
            },
            error = function(e) {
                err_tbl <- tibble(
                    surrogate_id = i,
                    effect_var = effect_var,
                    cause_var = cause_var,
                    edge_tag = tag,
                    error_message = conditionMessage(e),
                    saved_at = as.character(Sys.time())
                )
                surrogate_errors[[length(surrogate_errors) + 1L]] <<- err_tbl
                err_path <- file.path(surrogate_run_dir, paste0("surrogate_", surrogate_file_id(i), "_ERROR.csv"))
                save_csv(err_tbl, err_path)
                log_msg(
                    "Surrogate UIC failed for ", tag, " surrogate ", i,
                    "; error saved and continuing: ", conditionMessage(e),
                    level = "WARN"
                )
                log_progress("uic_surrogate_error", tag, paste0("surrogate ", i, ": ", conditionMessage(e)), err_path)
                NULL
            }
        )
        if (is.null(sres)) {
            next
        }
        surrogate_success <- surrogate_success + 1L
        surr_ete[, i] <- sres$ete[match(obs$tp, sres$tp)]
        if (isTRUE(config$save_each_surrogate)) {
            sres_out <- sres %>%
                mutate(
                    surrogate_id = i,
                    effect_var = effect_var,
                    cause_var = cause_var,
                    edge_tag = tag
                ) %>%
                relocate(surrogate_id, effect_var, cause_var, edge_tag)
            surrogate_result_path <- file.path(surrogate_run_dir, paste0("surrogate_", surrogate_file_id(i), "_uic.csv"))
            save_csv(sres_out, surrogate_result_path)
            log_progress(
                "uic_surrogate_saved",
                tag,
                paste0("surrogate ", i, "/", config$num_surr),
                surrogate_result_path
            )
        }
        if (isTRUE(config$save_intermediate) && (i %% config$checkpoint_every == 0 || i == config$num_surr)) {
            completed_matrix <- surr_ete[, seq_len(i), drop = FALSE]
            completed_tbl <- tibble(tp = obs$tp) %>%
                bind_cols(as_tibble(completed_matrix, .name_repair = "minimal"))
            progress_tbl <- tibble(
                effect_var = effect_var,
                cause_var = cause_var,
                edge_tag = tag,
                completed_surrogates = i,
                total_surrogates = config$num_surr,
                elapsed = format_elapsed(pair_start),
                saved_at = as.character(Sys.time())
            )
            save_csv(completed_tbl, file.path(pair_checkpoint_dir, paste0(tag, "_surrogate_ete_matrix_checkpoint_latest.csv")))
            save_csv(completed_tbl, file.path(pair_checkpoint_dir, paste0(tag, "_surrogate_ete_matrix_after_", surrogate_file_id(i), ".csv")))
            save_csv(progress_tbl, file.path(pair_checkpoint_dir, paste0(tag, "_surrogate_progress_latest.csv")))
            log_msg("Saved surrogate checkpoint for ", tag, ": ", i, "/", config$num_surr)
        }
        if (i %in% surr_progress) {
            log_msg(
                "Surrogate UIC progress for ", cause_var, " -> ", effect_var,
                ": ", i, "/", config$num_surr,
                " (", round(100 * i / config$num_surr), "%; elapsed ",
                format_elapsed(pair_start), ")"
            )
        }
    }

    log_msg("Calculating surrogate quantiles and empirical p-values for ", cause_var, " -> ", effect_var)
    if (length(surrogate_errors) > 0) {
        surrogate_error_path <- file.path(pair_checkpoint_dir, paste0(tag, "_surrogate_errors.csv"))
        save_csv(bind_rows(surrogate_errors), surrogate_error_path)
        log_progress(
            "uic_surrogate_errors_saved",
            tag,
            paste0(length(surrogate_errors), " failed surrogate run(s)"),
            surrogate_error_path
        )
    }
    if (surrogate_success < 1) {
        save_status(
            stage = "uic_surrogate",
            status = "failed",
            detail = paste0("All surrogate UIC runs failed for ", tag),
            path = pair_checkpoint_dir
        )
        stop("All surrogate UIC runs failed for ", tag, ". Inspect ", pair_checkpoint_dir, call. = FALSE)
    }

    qfun <- function(x, p) {
        x <- x[is.finite(x)]
        if (length(x) == 0) {
            return(NA_real_)
        }
        as.numeric(quantile(x, probs = p, names = FALSE))
    }
    q90 <- apply(surr_ete, 1, qfun, p = 0.90)
    q95 <- apply(surr_ete, 1, qfun, p = 0.95)
    q975 <- apply(surr_ete, 1, qfun, p = 0.975)
    q99 <- apply(surr_ete, 1, qfun, p = 0.99)

    p_emp <- vapply(seq_len(nrow(obs)), function(r) {
        valid <- surr_ete[r, is.finite(surr_ete[r, ])]
        if (length(valid) == 0 || !is.finite(obs$ete[r])) {
            return(NA_real_)
        }
        (1 + sum(valid >= obs$ete[r], na.rm = TRUE)) / (1 + length(valid))
    }, numeric(1))

    max_surr <- apply(surr_ete, 2, function(x) {
        x <- x[is.finite(x)]
        if (length(x) == 0) {
            return(NA_real_)
        }
        max(x)
    })
    valid_max_surr <- max_surr[is.finite(max_surr)]
    q95_global <- qfun(valid_max_surr, 0.95)
    p_global <- vapply(obs$ete, function(x) {
        if (length(valid_max_surr) == 0 || !is.finite(x)) {
            return(NA_real_)
        }
        (1 + sum(valid_max_surr >= x, na.rm = TRUE)) / (1 + length(valid_max_surr))
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
            n_surr_success = surrogate_success,
            n_surr_failed = length(surrogate_errors),
            sig_pointwise = .data$p_emp < config$alpha & .data$ete > .data$q95,
            sig_global = .data$p_global < config$alpha & .data$ete > .data$q95_global
        )
    log_msg(
        "UIC significance summary for ", cause_var, " -> ", effect_var,
        ": pointwise significant tp = ", sum(res$sig_pointwise, na.rm = TRUE),
        ", global significant tp = ", sum(res$sig_global, na.rm = TRUE)
    )

    uic_result_path <- file.path(paths$uic_tables, paste0(tag, "_uic_surrogate_corrected.csv"))
    uic_surr_matrix_path <- file.path(paths$uic_tables, paste0(tag, "_surrogate_ete_matrix.csv"))
    save_csv(res, uic_result_path)
    save_csv(as.data.frame(surr_ete), uic_surr_matrix_path)
    log_msg("Saved UIC tables for ", tag)
    log_progress(
        "uic_pair_completed",
        tag,
        paste0("successful surrogates: ", surrogate_success, "/", config$num_surr),
        uic_result_path
    )

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
    log_msg(
        "UIC completed: ", cause_var, " -> ", effect_var,
        " (elapsed ", format_elapsed(pair_start), ")"
    )

    res
}

uic_pair_grid <- tidyr::crossing(effect_var = vars, cause_var = vars) %>%
    filter(.data$effect_var != .data$cause_var)
log_msg("UIC pair plan: ", nrow(uic_pair_grid), " directed subtype pairs")

uic_results_list <- vector("list", nrow(uic_pair_grid))
for (i in seq_len(nrow(uic_pair_grid))) {
    uic_results_list[[i]] <- run_uic_pair(
        df_model = df_model,
        effect_var = uic_pair_grid$effect_var[[i]],
        cause_var = uic_pair_grid$cause_var[[i]],
        config = config,
        paths = paths,
        pair_id = i,
        total_pairs = nrow(uic_pair_grid)
    )
    if (isTRUE(config$save_intermediate)) {
        uic_partial <- bind_rows(uic_results_list[seq_len(i)])
        uic_partial_latest_path <- file.path(paths$uic_checkpoints, "uic_all_pairs_raw_partial_latest.csv")
        save_csv(uic_partial, uic_partial_latest_path)
        save_csv(
            uic_partial,
            file.path(paths$uic_checkpoints, paste0("uic_all_pairs_raw_partial_after_pair_", sprintf("%02d", i), ".csv"))
        )
        log_msg("Saved cumulative raw UIC checkpoint after pair ", i, "/", nrow(uic_pair_grid))
        log_progress(
            "uic_all_pairs_partial_saved",
            tag = paste0(uic_pair_grid$cause_var[[i]], "_to_", uic_pair_grid$effect_var[[i]]),
            detail = paste0("completed ", i, "/", nrow(uic_pair_grid), " UIC pairs"),
            path = uic_partial_latest_path
        )
    }
}
uic_all <- bind_rows(uic_results_list)
rm(uic_results_list)

log_msg("Skipping FDR correction across all UIC lag results; using unadjusted global empirical surrogate p-values")
uic_all <- uic_all %>%
    mutate(
        # No Benjamini-Hochberg/FDR correction is applied in this variant.
        # Primary UIC significance is the pair-wise global max-statistic surrogate test:
        #   p_global < alpha AND observed ETE > q95_global.
        sig_primary = .data$p_global < config$alpha & .data$ete > .data$q95_global
    )
if (isTRUE(config$save_intermediate)) {
    save_csv(uic_all, file.path(paths$uic_checkpoints, "uic_all_pairs_noFDR_checkpoint.csv"))
    log_msg("Saved no-FDR all-pair UIC checkpoint")
}

save_csv(uic_all, file.path(paths$tables, "uic_all_pairs_surrogate_corrected_noFDR.csv"))
log_msg("Saved all-pair no-FDR UIC table: ", file.path(paths$tables, "uic_all_pairs_surrogate_corrected_noFDR.csv"))

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
        pval = .data$p_global,
        te = ifelse(is.na(.data$te), .data$ete, .data$te)
    )

selected_links_path <- file.path(paths$tables, "uic_selected_links_for_MDR_noFDR.csv")
save_csv(selected_links, selected_links_path)
if (isTRUE(config$save_intermediate)) {
    save_csv(selected_links, file.path(paths$uic_checkpoints, "uic_selected_links_for_MDR_noFDR_checkpoint.csv"))
}
log_msg("Selected UIC links entering MDR: ", nrow(selected_links))
log_progress(
    "uic_selected_links_saved",
    tag = "noFDR",
    detail = paste0(nrow(selected_links), " selected link(s)"),
    path = selected_links_path
)
if (nrow(selected_links) > 0) {
    print(selected_links %>% select(effect_var, cause_var, tp, lag_weeks, ete, p_global))
}

analysis_has_mdr_links <- nrow(selected_links) > 0

if (!analysis_has_mdr_links) {
    no_link_msg <- paste0(
        "No UIC links passed the no-FDR seasonal-surrogate max-statistic criterion ",
        "(p_global < ", config$alpha, " and ete > q95_global). MDR S-map and manuscript Figure 4/5 outputs were skipped."
    )
    log_msg(no_link_msg, level = "WARN")

    no_link_status_path <- file.path(paths$tables, "analysis_completed_no_MDR_links.csv")
    save_csv(
        tibble(
            status = "completed_no_mdr_links",
            reason = no_link_msg,
            uic_all_pairs_table = file.path(paths$tables, "uic_all_pairs_surrogate_corrected_noFDR.csv"),
            selected_links_table = selected_links_path,
            saved_at = as.character(Sys.time())
        ),
        no_link_status_path
    )
    save_status("uic_selection", "completed_no_mdr_links", no_link_msg, no_link_status_path)

    log_msg("Building diagnostic UIC figure for all directed pairs because no links were selected")
    uic_plot_data_all <- uic_all %>%
        mutate(edge = edge_label(.data$cause_var, .data$effect_var))
    fig3_all <- ggplot(uic_plot_data_all, aes(x = .data$tp)) +
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
            title = "Diagnostic UIC plot: no no-FDR links passed the primary criterion"
        ) +
        theme_pub(13)
    fig3_all_path <- file.path(paths$fig, "Figure3_UIC_all_pairs_no_selected_links_noFDR.tiff")
    ggsave(fig3_all_path, fig3_all,
        width = 13, height = 9, dpi = config$dpi, compression = "lzw"
    )
    log_msg("Saved diagnostic UIC figure: ", fig3_all_path)
    log_progress("uic_diagnostic_figure_saved", "noFDR_no_links", no_link_msg, fig3_all_path)

    no_link_text <- c(
        "No MDR S-map links were selected in this run.",
        "",
        no_link_msg,
        "",
        paste0("Inspect the full UIC table: ", file.path(paths$tables, "uic_all_pairs_surrogate_corrected_noFDR.csv")),
        paste0("Inspect the diagnostic UIC figure: ", fig3_all_path)
    )
    writeLines(no_link_text, file.path(paths$manuscript, "NO_MDR_LINKS_noFDR.txt"))
}

if (analysis_has_mdr_links) {
    # Combined UIC figure for manuscript Figure 3.
    log_msg("Building combined UIC Figure 3 from selected links")
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
            title = "Figure 3. UIC with 52-week seasonal surrogate thresholds (no FDR correction)"
        ) +
        theme_pub(13)

    ggsave(file.path(paths$fig, "Figure3_UIC_seasonal_surrogate_selected_links_noFDR.tiff"), fig3,
        width = 13, height = 7, dpi = config$dpi, compression = "lzw"
    )
    log_msg("Saved Figure 3: ", file.path(paths$fig, "Figure3_UIC_seasonal_surrogate_selected_links_noFDR.tiff"))

    # -----------------------------
    # 6. MDR S-map per effect variable
    # -----------------------------
    log_section("6. MDR S-map per effect variable")

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

    run_mdr_for_effect <- function(effect_var, selected_links, df_model, df_log, config, paths, effect_id = NULL, total_effects = NULL) {
        effect_start <- Sys.time()
        effect_prefix <- if (!is.null(effect_id) && !is.null(total_effects)) {
            paste0("effect ", effect_id, "/", total_effects, ": ")
        } else {
            ""
        }
        log_msg("MDR S-map started: ", effect_prefix, "effect = ", effect_var)
        log_progress("mdr_effect_started", effect_var, paste0(effect_prefix, "effect = ", effect_var))

        links_eff <- selected_links %>% filter(.data$effect_var == effect_var)
        if (nrow(links_eff) == 0) {
            log_msg("No selected UIC links for effect = ", effect_var, "; skipping MDR", level = "WARN")
            log_progress("mdr_effect_skipped", effect_var, "No selected UIC links")
            return(NULL)
        }
        log_msg("MDR links for ", effect_var, ": ", nrow(links_eff))
        print(links_eff %>% select(effect_var, cause_var, tp, lag_weeks, ete, p_global))

        # Estimate effect-variable embedding dimension for the effect-history part of block_mvd.
        log_msg("Running simplex to choose effect-history embedding for ", effect_var)
        simp_x <- rUIC::simplex(
            df_model,
            lib_var = effect_var,
            E = config$E_range,
            tp = config$smap_tp,
            num_surr = config$uic_internal_num_surr
        )
        simp_x <- as_tibble(simp_x)
        Ex <- simp_x[which.min(simp_x$rmse), "E", drop = TRUE]
        Ex <- max(1, as.integer(Ex))
        save_csv(simp_x, file.path(paths$mdr_tables, paste0(effect_var, "_simplex_for_effect_embedding.csv")))
        log_msg("Simplex completed for ", effect_var, ": selected Ex = ", Ex)

        # Critical: make_block_mvd() internally filters on the column named 'pval'.
        # Therefore pval is replaced by the unadjusted global surrogate p-value for MDR block construction.
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
                pval = .data$p_global,
                q95_global = .data$q95_global,
                lag_weeks = .data$lag_weeks,
                edge = .data$edge
            )

        save_csv(uic_for_block, file.path(paths$mdr_tables, paste0(effect_var, "_uic_links_entering_MDR.csv")))
        log_msg("Saved MDR input UIC links for ", effect_var)

        # tp_adjust keeps the UIC lag interpretation consistent with the one-step S-map.
        # Example: UIC tp = -9 and S-map tp = 1 -> use cause_tp-8, so the target effect(t+1)
        # is paired with cause(t-8), preserving a 9-week cause-to-effect lag.
        log_msg("Building MDR block_mvd for ", effect_var, " using make_block_mvd()")
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

        block_mvd_path <- file.path(paths$mdr_tables, paste0(effect_var, "_block_mvd.csv"))
        save_csv(as_tibble(block_mvd), block_mvd_path)
        log_msg(
            "MDR block_mvd created for ", effect_var,
            ": ", nrow(block_mvd), " rows x ", ncol(block_mvd), " columns"
        )
        log_progress("mdr_block_saved", effect_var, paste0(ncol(block_mvd), " block columns"), block_mvd_path)

        mvd_E <- min(config$mdr_E, ncol(block_mvd))
        if (mvd_E < 1) stop("mvd_E became < 1; inspect block_mvd.", call. = FALSE)
        log_msg(
            "Computing multiview distance for ", effect_var,
            " (E = ", mvd_E, ", n_ssr = ", config$n_ssr, ", k = ", config$k, ")"
        )

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
        log_msg(
            "Multiview distance completed for ", effect_var,
            ": matrix ", nrow(multiview_dist), " x ", ncol(multiview_dist)
        )

        param_grid <- expand.grid(
            theta = config$theta_grid,
            lambda = config$lambda_grid,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
        log_msg(
            "Starting MDR theta/lambda search for ", effect_var,
            ": ", nrow(param_grid), " parameter combinations"
        )

        effect_checkpoint_dir <- safe_dir(file.path(paths$mdr_checkpoints, effect_var))
        param_results <- vector("list", nrow(param_grid))
        for (param_id in seq_len(nrow(param_grid))) {
            theta <- param_grid$theta[[param_id]]
            lambda <- param_grid$lambda[[param_id]]
            log_msg(
                "MDR parameter search ", effect_var,
                ": ", param_id, "/", nrow(param_grid),
                " theta = ", theta, ", lambda = ", lambda
            )
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
            fit_stats <- as_tibble(fit$stats) %>% mutate(theta = theta, lambda = lambda)
            param_results[[param_id]] <- fit_stats
            log_msg(
                "MDR parameter result ", effect_var,
                ": theta = ", theta,
                ", rho = ", signif(fit_stats$rho[[1]], 4),
                ", rmse = ", signif(fit_stats$rmse[[1]], 4)
            )
            if (isTRUE(config$save_intermediate)) {
                param_checkpoint_path <- file.path(effect_checkpoint_dir, paste0(effect_var, "_MDR_parameter_", sprintf("%03d", param_id), ".csv"))
                save_csv(fit_stats, param_checkpoint_path)
                param_partial <- bind_rows(param_results[seq_len(param_id)]) %>%
                    select(any_of(c("N", "theta", "lambda", "rho", "mae", "rmse")), everything())
                save_csv(param_partial, file.path(effect_checkpoint_dir, paste0(effect_var, "_MDR_parameter_search_partial_latest.csv")))
                log_msg("Saved MDR parameter-search checkpoint for ", effect_var, ": ", param_id, "/", nrow(param_grid))
                log_progress(
                    "mdr_parameter_checkpoint_saved",
                    effect_var,
                    paste0("parameter ", param_id, "/", nrow(param_grid)),
                    param_checkpoint_path
                )
            }
        }
        param_res <- bind_rows(param_results) %>%
            select(any_of(c("N", "theta", "lambda", "rho", "mae", "rmse")), everything())
        rm(param_results)

        param_res_path <- file.path(paths$mdr_tables, paste0(effect_var, "_MDR_parameter_search.csv"))
        save_csv(param_res, param_res_path)
        log_msg("Saved MDR parameter search table for ", effect_var)
        log_progress("mdr_parameter_search_saved", effect_var, paste0(nrow(param_res), " parameter rows"), param_res_path)

        best <- param_res %>%
            filter(is.finite(.data$rmse)) %>%
            arrange(.data$rmse, desc(.data$rho)) %>%
            slice(1)
        if (nrow(best) == 0) stop("No valid MDR parameter result for ", effect_var, call. = FALSE)
        log_msg(
            "Best MDR parameters for ", effect_var,
            ": theta = ", best$theta[[1]],
            ", lambda = ", best$lambda[[1]],
            ", rho = ", signif(best$rho[[1]], 4),
            ", rmse = ", signif(best$rmse[[1]], 4)
        )

        log_msg("Running final MDR S-map with coefficient export for ", effect_var)
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
        log_msg("Saved final MDR fit RDS for ", effect_var)

        pred <- as_tibble(final_fit$model_output) %>%
            mutate(Date = pred_dates(., df_log), effect_var = effect_var) %>%
            relocate(Date, effect_var)
        coefs <- as_tibble(final_fit$smap_coefficients) %>%
            mutate(Date = coef_dates(., df_log), effect_var = effect_var) %>%
            relocate(Date, effect_var)

        save_csv(pred, file.path(paths$mdr_tables, paste0(effect_var, "_MDR_prediction.csv")))
        save_csv(coefs, file.path(paths$mdr_tables, paste0(effect_var, "_MDR_coefficients_raw.csv")))
        log_msg(
            "Saved MDR predictions and raw coefficients for ", effect_var,
            ": predictions = ", nrow(pred), " rows; coefficients = ", nrow(coefs), " rows"
        )

        cmap <- coef_name_map(block_mvd, coefs)
        save_csv(cmap, file.path(paths$mdr_tables, paste0(effect_var, "_MDR_coefficient_column_map.csv")))
        if (any(!cmap$exists)) {
            warning(
                "Some expected coefficient columns were not found for ", effect_var, ": ",
                paste(cmap$coef_col[!cmap$exists], collapse = ", ")
            )
        }
        log_msg("Mapped MDR coefficient columns for ", effect_var)

        log_msg("Extracting edge-specific MDR coefficients for ", effect_var)
        edge_coef <- purrr::map_dfr(seq_len(nrow(uic_for_block)), function(i) {
            link <- uic_for_block[i, ]
            adjusted_tp <- link$tp + config$smap_tp
            block_col <- paste0(link$cause_var, "_tp", adjusted_tp)
            coef_col <- cmap$coef_col[match(block_col, cmap$block_col)]

            if (length(coef_col) == 0 || is.na(coef_col) || !(coef_col %in% names(coefs))) {
                warning("Coefficient for ", link$edge, " not found. Expected block column: ", block_col)
                return(tibble())
            }

            coef_values <- coefs[[coef_col]]
            tibble(
                Date = coefs$Date,
                effect_var = link$effect_var,
                cause_var = link$cause_var,
                edge = link$edge,
                lag_weeks = link$lag_weeks,
                uic_tp = link$tp,
                mdr_block_col = block_col,
                coef_col = coef_col,
                coef_value = coef_values
            )
        })

        save_csv(edge_coef, file.path(paths$mdr_tables, paste0(effect_var, "_MDR_edge_coefficients_long.csv")))
        log_msg("Saved edge coefficient long table for ", effect_var, ": ", nrow(edge_coef), " rows")
        if (nrow(edge_coef) == 0) {
            stop(
                "No MDR edge coefficients were extracted for ", effect_var,
                ". Check the block_mvd column names and coefficient map table.",
                call. = FALSE
            )
        }

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
                mae = if ("mae" %in% names(best)) best$mae[[1]] else NA_real_,
                rmse = best$rmse[[1]],
                regularized = config$ridge_regularized_mdr
            )

        edge_summary_path <- file.path(paths$mdr_tables, paste0(effect_var, "_MDR_edge_summary.csv"))
        save_csv(edge_summary, edge_summary_path)
        log_msg("Saved edge summary for ", effect_var, ": ", nrow(edge_summary), " rows")
        log_progress("mdr_edge_summary_saved", effect_var, paste0(nrow(edge_summary), " edge summary rows"), edge_summary_path)
        print(edge_summary %>% select(effect_var, cause_var, lag_weeks, coef_median, prop_positive, polarity, theta, rho, rmse))

        # Per-effect diagnostic plot: observed vs predicted.
        log_msg("Saving MDR observed-vs-predicted plot for ", effect_var)
        p_obs_pred <- ggplot(pred, aes(x = .data$obs, y = .data$pred)) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
            geom_point(alpha = 0.65, size = 1.7, na.rm = TRUE) +
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
            log_msg("Saving MDR coefficient time-series plot for ", effect_var)
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
        log_msg("MDR S-map completed for ", effect_var, " (elapsed ", format_elapsed(effect_start), ")")
        log_progress("mdr_effect_completed", effect_var, paste0("elapsed ", format_elapsed(effect_start)))

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

    mdr_effect_vars <- unique(selected_links$effect_var)
    log_msg("MDR effect plan: ", length(mdr_effect_vars), " effect variable(s): ", paste(mdr_effect_vars, collapse = ", "))
    mdr_results <- purrr::map(seq_along(mdr_effect_vars), function(i) {
        effect_var <- mdr_effect_vars[[i]]
        run_mdr_for_effect(
            effect_var = effect_var,
            selected_links = selected_links,
            df_model = df_model,
            df_log = df_log,
            config = config,
            paths = paths,
            effect_id = i,
            total_effects = length(mdr_effect_vars)
        )
    })
    names(mdr_results) <- mdr_effect_vars
    mdr_results <- purrr::compact(mdr_results)
    if (length(mdr_results) == 0) {
        stop("No MDR S-map results were generated. Inspect selected UIC links and MDR logs.", call. = FALSE)
    }
    save_rds(mdr_results, file.path(paths$mdr_tables, "MDR_results_all_effects.rds"))
    log_msg("Saved all MDR results RDS: ", file.path(paths$mdr_tables, "MDR_results_all_effects.rds"))

    edge_coef_all <- purrr::map_dfr(mdr_results, "edge_coef")
    edge_summary_all <- purrr::map_dfr(mdr_results, "edge_summary")
    log_msg(
        "Combined MDR outputs: edge coefficients = ", nrow(edge_coef_all),
        " rows; edge summaries = ", nrow(edge_summary_all), " rows"
    )

    # Join UIC information onto MDR summary.
    log_msg("Joining UIC information onto MDR summary")
    edge_summary_all <- edge_summary_all %>%
        left_join(
            selected_links %>%
                select(
                    effect_var, cause_var, E, tp, ete,
                    q95_global, p_global
                ),
            by = c("effect_var", "cause_var")
        ) %>%
        rename(
            uic_E = E,
            uic_tp_selected = tp,
            uic_ete = ete,
            uic_q95_global = q95_global,
            uic_p_global = p_global
        ) %>%
        arrange(.data$effect_var, .data$cause_var)

    edge_coef_all_path <- file.path(paths$tables, "Figure5_MDR_coefficients_long.csv")
    edge_summary_all_path <- file.path(paths$tables, "Table_MDR_edge_summary.csv")
    save_csv(edge_coef_all, edge_coef_all_path)
    save_csv(edge_summary_all, edge_summary_all_path)
    log_msg("Saved combined MDR coefficient and summary tables")
    log_progress(
        "mdr_combined_tables_saved",
        "all_effects",
        paste0(nrow(edge_coef_all), " coefficient rows; ", nrow(edge_summary_all), " summary rows"),
        edge_summary_all_path
    )

    # Compatibility export for downstream figure scripts. This does not overwrite old result folders.
    log_msg("Creating wide MDR coefficient export for downstream Figure 5 scripts")
    coef_wide <- edge_coef_all %>%
        mutate(edge_col = paste0(.data$cause_var, "_cause_", .data$effect_var)) %>%
        group_by(.data$edge_col) %>%
        mutate(row_id = row_number()) %>%
        ungroup() %>%
        select(row_id, edge_col, coef_value) %>%
        pivot_wider(names_from = edge_col, values_from = coef_value) %>%
        select(-row_id)
    coef_wide_path <- file.path(paths$tables, "MDR_coefficients_all_wide_for_Figure5.csv")
    save_csv(coef_wide, coef_wide_path)
    log_msg("Saved wide MDR coefficient table: ", coef_wide_path)
    log_progress("mdr_wide_coefficients_saved", "all_effects", paste0(ncol(coef_wide), " edge columns"), coef_wide_path)

    # -----------------------------
    # 7. Manuscript-ready figures
    # -----------------------------
    log_section("7. Manuscript-ready figures")

    log_msg("Building Figure 5 coefficient distribution plot")
    fig5 <- edge_coef_all %>%
        mutate(edge = forcats::fct_reorder(.data$edge, .data$coef_value, .fun = median, na.rm = TRUE, .na_rm = TRUE)) %>%
        ggplot(aes(x = .data$edge, y = .data$coef_value)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_violin(alpha = 0.45, trim = FALSE, na.rm = TRUE) +
        geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.75, na.rm = TRUE) +
        ggforce::geom_sina(alpha = 0.35, size = 0.8, na.rm = TRUE) +
        coord_flip() +
        labs(
            title = "Figure 5. MDR S-map coefficients for no-FDR UIC-selected subtype interactions",
            x = NULL,
            y = "MDR S-map coefficient on log(x+1) scale"
        ) +
        theme_pub(14)

    ggsave(file.path(paths$fig, "Figure5_MDR_Smap_coefficients.tiff"), fig5,
        width = 10, height = max(5, 1.5 + length(unique(edge_coef_all$edge)) * 1.0),
        dpi = config$dpi, compression = "lzw"
    )
    log_msg("Saved Figure 5: ", file.path(paths$fig, "Figure5_MDR_Smap_coefficients.tiff"))

    # Figure 4 network summary.
    log_msg("Building Figure 4 network summary")
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
    log_msg("Saved Figure 4 edge data: ", file.path(paths$tables, "Figure4_network_edges_MDR.csv"))

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
            title = "Figure 4. Directional subtype interaction network from no-FDR UIC and MDR S-map",
            linetype = "median coefficient"
        ) +
        theme_void(base_size = 14) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")

    ggsave(file.path(paths$fig, "Figure4_UIC_MDR_network.tiff"), fig4,
        width = 8, height = 6, dpi = config$dpi, compression = "lzw"
    )
    log_msg("Saved Figure 4: ", file.path(paths$fig, "Figure4_UIC_MDR_network.tiff"))

    # -----------------------------
    # 8. Manuscript-ready tables and text snippets
    # -----------------------------
    log_section("8. Manuscript-ready tables and text snippets")

    log_msg("Building manuscript-ready main results table")
    table_main <- edge_summary_all %>%
        transmute(
            Cause = subtype_label(.data$cause_var),
            Effect = subtype_label(.data$effect_var),
            `UIC lag, weeks` = .data$lag_weeks,
            `UIC E` = .data$uic_E,
            `UIC ETE` = round(.data$uic_ete, 4),
            `UIC global p` = signif(.data$uic_p_global, 3),
            `MDR theta` = signif(.data$theta, 3),
            `MDR rho` = signif(.data$rho, 3),
            `MDR RMSE` = signif(.data$rmse, 3),
            `Mean MDR coef` = round(.data$coef_mean, 4),
            `Median MDR coef` = round(.data$coef_median, 4),
            `IQR` = paste0(round(.data$coef_q25, 4), " to ", round(.data$coef_q75, 4)),
            `Positive proportion` = round(.data$prop_positive, 3),
            Polarity = .data$polarity
        )

    table_main_path <- file.path(paths$tables, "Table1_UIC_MDR_main_results_noFDR.csv")
    save_csv(table_main, table_main_path)
    log_msg("Saved Table 1 main results: ", table_main_path)
    log_progress("manuscript_table_saved", "Table1", paste0(nrow(table_main), " rows"), table_main_path)
    print(table_main)

    methods_text <- c(
        "Suggested Methods replacement for the S-map paragraph:",
        "",
        paste0(
            "Interaction strength and polarity were estimated using multiview-distance S-map (MDR S-map). ",
            "For each effect subtype, UIC was first performed across the remaining subtypes over lags of ",
            paste(range(config$tp_range), collapse = " to "),
            " using rUIC::uic.optimal(), which selects an embedding dimension before returning lag-wise UIC values. ",
            "UIC links were retained for MDR S-map only when the observed effective transfer entropy exceeded ",
            "the 95th percentile of the 52-week seasonal surrogate max-statistic distribution and the empirical ",
            "global empirical surrogate p-value was below the prespecified alpha threshold without Benjamini-Hochberg/FDR correction. ",
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
        "Suggested Results table is exported as tables/Table1_UIC_MDR_main_results_noFDR.csv. Replace numerical values in the manuscript after running the script on the final dataset."
    )
    writeLines(methods_text, file.path(paths$manuscript, "MDR_methods_results_figure_legend_text.txt"))
    log_msg("Saved manuscript text snippets: ", file.path(paths$manuscript, "MDR_methods_results_figure_legend_text.txt"))
}

# -----------------------------
# 9. Session information
# -----------------------------
log_section("9. Session information")

log_msg("Writing sessionInfo.txt")
final_status <- if (isTRUE(analysis_has_mdr_links)) {
    "completed_with_mdr"
} else {
    "completed_no_mdr_links"
}
write_session_info(status = final_status)
log_msg("Saved session info: ", file.path(paths$root, "sessionInfo.txt"))
save_status(
    stage = "analysis",
    status = final_status,
    detail = "Pipeline finished. See progress/progress_log.csv for incremental outputs.",
    path = paths$root
)

log_msg(
    "Done. Results written to: ",
    normalizePath(paths$root, winslash = "/", mustWork = FALSE),
    " (total elapsed ", format_elapsed(analysis_start_time), ")"
)
