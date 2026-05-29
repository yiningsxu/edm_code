# Resume only the MDR S-map part of flu_subtype_mdrsmap.R.
#
# Default behavior:
#   - Finds the latest result/FluSub_JP/*_UICoptimal_MDR_noFDR folder
#     containing the no-FDR UIC selected-link table.
#   - Reuses prepared_log1p_timeseries.csv and uic_selected_links_for_MDR_noFDR.csv.
#   - Writes MDR outputs to a new sibling folder named *_MDR_resume_YYYYMMDD_HHMMSS.
#
# Optional environment variables:
#   FLU_SUBTYPE_RESUME_DIR : source result folder with completed UIC outputs.
#   FLU_SUBTYPE_OUT_DIR    : destination folder for resumed MDR outputs.
#   FLU_SUBTYPE_N_SSR, FLU_SUBTYPE_THETA_GRID, etc. override values from
#                            the source run_config.csv.

rm(list = ls())

current_script_path <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- args[startsWith(args, "--file=")][1]
    if (!is.na(file_arg)) {
        file_arg <- sub("^--file=", "", file_arg)
        if (file.exists(file_arg)) {
            return(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
        }
    }

    for (frame in rev(sys.frames())) {
        ofile <- frame$ofile
        if (!is.null(ofile) && nzchar(ofile) && file.exists(ofile)) {
            return(normalizePath(ofile, winslash = "/", mustWork = TRUE))
        }
    }

    NA_character_
}

find_edm_code_root <- function(start_dir) {
    current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
    repeat {
        if (file.exists(file.path(current, ".edm_code_root"))) {
            return(current)
        }
        parent <- dirname(current)
        if (identical(parent, current)) {
            stop("Could not find edm_code root from: ", start_dir, call. = FALSE)
        }
        current <- parent
    }
}

find_latest_resume_dir <- function(workspace_dir) {
    candidates <- list.dirs(file.path(workspace_dir, "result", "FluSub_JP"),
        full.names = TRUE, recursive = FALSE
    )
    candidates <- candidates[grepl("_UICoptimal_MDR_noFDR$", basename(candidates))]
    selected_paths <- file.path(candidates, "tables", "uic_selected_links_for_MDR_noFDR.csv")
    prepared_paths <- file.path(candidates, "tables", "prepared_log1p_timeseries.csv")
    candidates <- candidates[file.exists(selected_paths) & file.exists(prepared_paths)]

    if (length(candidates) == 0) {
        stop(
            "No completed no-FDR UIC result folder found. Set FLU_SUBTYPE_RESUME_DIR.",
            call. = FALSE
        )
    }

    selected_paths <- file.path(candidates, "tables", "uic_selected_links_for_MDR_noFDR.csv")
    selected_rows <- vapply(selected_paths, function(path) {
        max(0L, length(readLines(path, warn = FALSE)) - 1L)
    }, integer(1))
    candidates_with_links <- candidates[selected_rows > 0]

    if (length(candidates_with_links) == 0) {
        stop(
            "No no-FDR UIC result folder has selected MDR links. Candidate row counts: ",
            paste0(basename(candidates), "=", selected_rows, collapse = ", "),
            ". Set FLU_SUBTYPE_RESUME_DIR to a folder whose uic_selected_links_for_MDR_noFDR.csv has at least one data row.",
            call. = FALSE
        )
    }

    candidates <- candidates_with_links
    candidates[which.max(file.info(candidates)$mtime)]
}

config_value <- function(config_tbl, name) {
    value <- config_tbl$value[match(name, config_tbl$name)]
    if (length(value) == 0 || is.na(value)) {
        return(NA_character_)
    }
    value
}

set_env_if_unset <- function(env_name, value) {
    current <- Sys.getenv(env_name, unset = NA_character_)
    if (!is.na(current) && nzchar(current)) {
        return(invisible(FALSE))
    }
    if (is.na(value) || !nzchar(value)) {
        return(invisible(FALSE))
    }
    do.call(Sys.setenv, setNames(list(value), env_name))
    invisible(TRUE)
}

script_path <- current_script_path()
script_dir <- if (!is.na(script_path)) dirname(script_path) else getwd()
edm_code_root <- find_edm_code_root(script_dir)
workspace_dir <- normalizePath(file.path(edm_code_root, ".."), winslash = "/", mustWork = TRUE)
setwd(workspace_dir)

primary_script <- file.path(edm_code_root, "projects", "jp_flu_subtypes", "scripts", "flu_subtype_mdrsmap.R")
if (!file.exists(primary_script)) {
    stop("Primary no-FDR script not found: ", primary_script, call. = FALSE)
}

resume_dir <- Sys.getenv("FLU_SUBTYPE_RESUME_DIR", unset = NA_character_)
if (is.na(resume_dir) || !nzchar(resume_dir)) {
    resume_dir <- find_latest_resume_dir(workspace_dir)
}
resume_dir <- normalizePath(resume_dir, winslash = "/", mustWork = TRUE)

run_config_path <- file.path(resume_dir, "tables", "run_config.csv")
if (file.exists(run_config_path)) {
    source_config <- read.csv(run_config_path, stringsAsFactors = FALSE, check.names = FALSE)
    env_map <- c(
        data_file = "FLU_SUBTYPE_DATA_FILE",
        E_range = "FLU_SUBTYPE_E_RANGE",
        tp_range = "FLU_SUBTYPE_TP_RANGE",
        tau = "FLU_SUBTYPE_TAU",
        alpha = "FLU_SUBTYPE_ALPHA",
        num_surr = "FLU_SUBTYPE_NUM_SURR",
        uic_internal_num_surr = "FLU_SUBTYPE_UIC_INTERNAL_NUM_SURR",
        season_period = "FLU_SUBTYPE_SEASON_PERIOD",
        random_seed = "FLU_SUBTYPE_RANDOM_SEED",
        save_intermediate = "FLU_SUBTYPE_SAVE_INTERMEDIATE",
        save_each_surrogate = "FLU_SUBTYPE_SAVE_EACH_SURROGATE",
        checkpoint_every = "FLU_SUBTYPE_CHECKPOINT_EVERY",
        smap_tp = "FLU_SUBTYPE_SMAP_TP",
        mdr_E = "FLU_SUBTYPE_MDR_E",
        n_ssr = "FLU_SUBTYPE_N_SSR",
        k = "FLU_SUBTYPE_K",
        theta_grid = "FLU_SUBTYPE_THETA_GRID",
        ridge_regularized_mdr = "FLU_SUBTYPE_RIDGE_REGULARIZED_MDR",
        lambda_grid = "FLU_SUBTYPE_LAMBDA_GRID",
        dpi = "FLU_SUBTYPE_DPI"
    )
    for (name in names(env_map)) {
        set_env_if_unset(env_map[[name]], config_value(source_config, name))
    }
}

out_dir <- Sys.getenv("FLU_SUBTYPE_OUT_DIR", unset = NA_character_)
if (is.na(out_dir) || !nzchar(out_dir)) {
    out_dir <- paste0(resume_dir, "_MDR_resume_", format(Sys.time(), "%Y%m%d_%H%M%S"))
}
Sys.setenv(FLU_SUBTYPE_OUT_DIR = out_dir)

source_lines <- readLines(primary_script, warn = FALSE)
header_end <- grep("^# 3\\. Data preparation", source_lines)[1] - 2L
section6_start <- grep("^    # 6\\. MDR S-map per effect variable", source_lines)[1] - 1L
section9_start <- grep("^# 9\\. Session information", source_lines)[1]
section6_end <- section9_start - 4L

if (anyNA(c(header_end, section6_start, section6_end)) || section6_start >= section6_end) {
    stop("Could not locate MDR section markers in: ", primary_script, call. = FALSE)
}

header_code <- paste(source_lines[seq_len(header_end)], collapse = "\n")
mdr_code <- paste(source_lines[section6_start:section6_end], collapse = "\n")

run_env <- new.env(parent = globalenv())
eval(parse(text = header_code), envir = run_env)

prepared_path <- file.path(resume_dir, "tables", "prepared_log1p_timeseries.csv")
selected_path <- file.path(resume_dir, "tables", "uic_selected_links_for_MDR_noFDR.csv")
uic_all_path <- file.path(resume_dir, "tables", "uic_all_pairs_surrogate_corrected_noFDR.csv")

if (!file.exists(prepared_path)) {
    stop("Prepared time-series table not found: ", prepared_path, call. = FALSE)
}
if (!file.exists(selected_path)) {
    stop("Selected UIC links table not found: ", selected_path, call. = FALSE)
}

df_log <- read.csv(prepared_path, stringsAsFactors = FALSE, check.names = FALSE)
df_log$Date <- as.Date(df_log$Date)
selected_links <- read.csv(selected_path, stringsAsFactors = FALSE, check.names = FALSE)

required_selected <- c(
    "effect_var", "cause_var", "E", "tp", "te", "ete", "pval_raw_uic",
    "p_global", "q95_global", "lag_weeks", "edge"
)
missing_selected <- setdiff(required_selected, names(selected_links))
if (length(missing_selected) > 0) {
    stop(
        "Selected UIC links table is missing required column(s): ",
        paste(missing_selected, collapse = ", "),
        call. = FALSE
    )
}
if (nrow(selected_links) == 0) {
    stop("Selected UIC links table has 0 rows; there is no MDR section to resume.", call. = FALSE)
}

missing_model_cols <- setdiff(run_env$config$subtype_vars, names(df_log))
if (length(missing_model_cols) > 0) {
    stop(
        "Prepared time-series table is missing subtype column(s): ",
        paste(missing_model_cols, collapse = ", "),
        call. = FALSE
    )
}

run_env$resume_dir <- resume_dir
run_env$df_log <- df_log
run_env$df_model <- df_log[, run_env$config$subtype_vars, drop = FALSE]
run_env$selected_links <- selected_links
run_env$analysis_has_mdr_links <- TRUE

run_env$log_section("3. MDR resume inputs")
run_env$log_msg("Resume source folder: ", resume_dir)
run_env$log_msg("Resume output folder: ", normalizePath(run_env$paths$root, winslash = "/", mustWork = FALSE))
run_env$log_msg("Loaded prepared time-series rows: ", nrow(df_log))
run_env$log_msg("Loaded selected UIC links: ", nrow(selected_links))

run_env$save_csv(df_log, file.path(run_env$paths$tables, "prepared_log1p_timeseries.csv"))
run_env$save_csv(selected_links, file.path(run_env$paths$tables, "uic_selected_links_for_MDR_noFDR.csv"))
if (file.exists(uic_all_path)) {
    uic_all <- read.csv(uic_all_path, stringsAsFactors = FALSE, check.names = FALSE)
    run_env$uic_all <- uic_all
    run_env$save_csv(uic_all, file.path(run_env$paths$tables, "uic_all_pairs_surrogate_corrected_noFDR.csv"))
}
run_env$save_status(
    stage = "mdr_resume_inputs",
    status = "loaded",
    detail = paste0("Loaded ", nrow(selected_links), " selected UIC link(s) from prior UIC run."),
    path = resume_dir
)

eval(parse(text = mdr_code), envir = run_env)

eval(parse(text = '
log_section("9. Session information")
log_msg("Writing sessionInfo.txt")
write_session_info(
    status = "completed_mdr_resume",
    notes = c(paste0("MDR resume source folder: ", resume_dir))
)
log_msg("Saved session info: ", file.path(paths$root, "sessionInfo.txt"))
save_status(
    stage = "mdr_resume",
    status = "completed_mdr_resume",
    detail = "MDR-only resume finished. UIC outputs were reused from the source folder.",
    path = paths$root
)
log_msg(
    "Done. MDR resume results written to: ",
    normalizePath(paths$root, winslash = "/", mustWork = FALSE),
    " (total elapsed ", format_elapsed(analysis_start_time), ")"
)
'), envir = run_env)
