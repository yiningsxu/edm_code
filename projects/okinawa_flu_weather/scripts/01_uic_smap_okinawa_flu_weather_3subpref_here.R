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
# install.packages("devtools")
# remotes::install_github("ha0ye/rEDM")
# remotes::install_github("yutakaos/rUIC")
# remotes::install_github("ong8181/macam")

pacman::p_load(
  lubridate, tidyverse, ISOweek, ggplot2, ggforce, glue, here,
  gridExtra, cowplot, rlang, rEDM, macamts, rUIC, sinaplot, ggExtra, ggdensity
)
# Load library
packageVersion("rEDM") # v 0.7.5
packageVersion("macam") # v 0.1.4 → v 0.1.5
# packageVersion("macamts") # v 0.1.4
packageVersion("rUIC") # v 0.9.12
# library(rEDM)
# vignette("rEDM-tutorial")
theme_set(theme_cowplot())

# ----------------------------------------------- #
## Set path
base_dir <- workspace_root()
# setwd handled by edm_code bootstrap
source(project_file("projects", "okinawa_flu_weather", "scripts", "functions_oknwFlu.R"))

base_path <- file.path("result", "oknw_Flu", as.character(Sys.Date()))
dir.create(base_path, recursive = TRUE, showWarnings = FALSE)

# # Analysis Parameters
# REGIONS <- c("Yaeyama", "mainIsland", "Miyako")
# WEATHER_TYPES <- c("ah","temp","humid","rain") # 解析したい気象変数を指定 (例: c("temp", "ah"))
# NUM_SURROGATES <- 2000
# EFFECT_VAR <- "flu"

## ------------------------------------------------------------------------------------------------------------ ##
# 2. 関数の定義 (Function Definitions)
## ------------------------------------------------------------------------------------------------------------ ##
make_dirs_for_region <- function(base_path, area) {
  dirs <- list(
    area = file.path(base_path, area),
    fig = file.path(base_path, area, "figure"),
    smap_fig = file.path(base_path, area, "figure", "smapFig"),
    raw_ts = file.path(base_path, area, "figure", "raw_ts"),
    uic_fig = file.path(base_path, area, "figure", "uicFig", "uic"),
    uic_surr = file.path(base_path, area, "figure", "uicFig", "uic_surr"),
    uic_fig_dir = file.path(base_path, area, "figure", "uicFig"),
    xmap = file.path(base_path, area, "figure", "crossMapping"),
    table = file.path(base_path, area, "table"),
    smap_table = file.path(base_path, area, "table", "smap"),
    smap_block = file.path(base_path, area, "table", "smap", "block"),
    smap_coef = file.path(base_path, area, "table", "smap", "coef"),
    smap_pred = file.path(base_path, area, "table", "smap", "pred_res"),
    uic_tbl = file.path(base_path, area, "table", "uic", "result"),
    uic_surr_tbl = file.path(base_path, area, "table", "uic", "surrogate_dt"),
    uic_tbl_dir = file.path(base_path, area, "table", "uic")
  )
  for (d in unique(unlist(dirs))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  return(dirs)
}

# # absolute humidity calculation (T in degC, RH in %)
# calc_absolute_humidity <- function(temp_c, rh_percent) {
#   es <- 6.1078 * 10^(7.5 * temp_c / (237.3 + temp_c)) # hPa
#   e <- (rh_percent / 100) * es
#   ah <- 216.7 * e / (temp_c + 273.15) # g/m^3
#   return(ah)
# }

# --- main processing function for a region ---
UIC_func <- function(data_input, region) {
  flu_col <- glue("flu_{region}")
  temp_col <- glue("temp_{region}")
  humid_col <- glue("humid_{region}")
  rain_col <- glue("rainfall_{region}")
  print(paste(flu_col, temp_col, humid_col, rain_col))

  # # 飽和水蒸気圧 [hPa]
  # es <- 6.1078 * 10^(7.5 * dt0[[temp_col]] / (237.3 + dt0[[temp_col]]))
  # # 実際の水蒸気圧 [hPa]
  # e <- (dt0[[humid_col]] / 100) * es
  # # 絶対湿度 [g/m^3]
  # dt0$ah_col <- 216.7 * e / (dt0[[temp_col]] + 273.15)

  # log(): base "e"
  df_log <- data_input %>%
    # 必要な列を選択し、分かりやすい名前に変更
    select(
      date,
      flu = all_of(flu_col),
      temp = all_of(temp_col),
      humid = all_of(humid_col),
      rain = all_of(rain_col)
    ) %>%
    # 絶対湿度(ah)を計算
    mutate(
      es = 6.1078 * 10^(7.5 * temp / (237.3 + temp)), # 飽和水蒸気圧
      e = (humid / 100) * es, # 実際の水蒸気圧
      ah = 216.7 * e / (temp + 273.15), # 絶対湿度
    ) %>%
    # 各変数を対数変換 (log(x+1))
    mutate(across(c(flu, temp, humid, ah, rain), ~ log(.x + 1)))
  head(df_log)

  # assign(paste0("df_flu", region, "_norm"),
  #        data.frame(Date = df0_flu$date,
  #                   flu = df0_flu[[flu_col]] / mean(df0_flu[[flu_col]], na.rm = TRUE),
  #                   temp = df0_flu[[temp_col]] / mean(df0_flu[[temp_col]], na.rm = TRUE),
  #                   humid = df0_flu[[humid_col]] / mean(df0_flu[[humid_col]], na.rm = TRUE),
  #                   rain = df0_flu[[rain_col]] / mean(df0_flu[[rain_col]], na.rm = TRUE))
  # )

  # uic_uniRegSmap(get(paste0("df_flu", region, "_log")), get(paste0("df_flu", region, "_norm")), region)
  # uic_uniRegSmap(paste0("df_flu", region, "_log"), paste0("df_flu", region, "_norm"), region)

  area <- region
  print(paste0("------------------------------ ", area, " ------------------------------"))

  # Create directories using the helper function
  dirs <- make_dirs_for_region(base_path, area)

  # Assign paths to variables used in the rest of the script
  # Ensure trailing slashes are added where the original code had them,
  # to maintain compatibility with paste0() calls in other functions.
  res_save_path <- paste0(dirs$area, "/")
  fig_date_path <- paste0(dirs$fig, "/")
  uic_plot_path <- paste0(dirs$uic_fig_dir, "/")
  uic_surr_tbl_path <- paste0(dirs$uic_tbl_dir, "/")
  xmap_plot_path <- paste0(dirs$xmap, "/")

  print(res_save_path)
  print(fig_date_path)
  print(uic_surr_tbl_path)

  # plot with the log ts
  df_log_long <- df_log %>%
    pivot_longer(cols = c(flu, temp, humid, ah, rain), names_to = "variable", values_to = "value")

  ggplot(df_log_long, aes(x = date, y = value, color = variable)) +
    geom_line() +
    geom_point() +
    labs(color = "Variable") +
    theme(text = element_text(size = 20)) +
    xlab("Date") +
    ylab("Log Transformed Incidence") +
    facet_wrap(~variable, scales = "free_y", ncol = 1)

  ggsave(file.path(base_path, area, "figure", "raw_ts", "log_ts.tiff"),
    units = "in", width = 15, height = 10, dpi = 300, compression = "lzw"
  )
  print("plot and save log time series figure --- done")

  ## -------------------------- Seasonal surrogate (UIC) -------------------------- ##
  # List of influenza subtypes
  weathertypes <- c("ah", "temp", "humid", "rain")
  # weathertypes <- c("ah")
  # Perform UIC analysis for all pairs of subtypes
  effect_var <- "flu"
  for (cause_var in weathertypes) {
    numSurr <- 2000
    # numSurr <- 5
    Surr_UIC(df_log, effect_var, cause_var, numSurr, uic_surr_tbl_path, uic_plot_path)
  }
}

## ------------------------------------------------------------------------------------------------------------ ##
# 3. 実行 (Execution)
## ------------------------------------------------------------------------------------------------------------ ##
# --- load data ---
dt0 <- read.csv(file.path("data", "flu_subpref_oknw", "oknw_flu_weather_241226.csv"), stringsAsFactors = FALSE)
# construct date from ISO week (assumes columns 'year' and 'week' exist)
dt0 <- dt0 %>%
  mutate(date = as.Date(ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1"))))

# --- run for multiple regions ---
regions <- c("mainIsland", "Miyako", "Yaeyama")
for (r in regions) {
  UIC_func(dt0, r)
}
