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

rm(list = ls())

pacman::p_load(
  lubridate,
  tidyverse,
  ISOweek,
  zoo,
  # rEDM,
  ggplot2,
  # ggforce,
  # glue,
  # stats,
  dplyr,
  # gridExtra,
  # cowplot,
  # rlang,
  # # macam,
  # macamts,
  # rUIC,
  # sinaplot,
  # ggExtra,
  # ggdensity,
  # paletteer,
  viridis
)

theme_set(theme_cowplot())


dt_ID <- read.csv(
  "data/oknw/oknwID_10to19_mainIsland_2islands_250314.csv"
)
dt_ID <- dt_ID %>%
  mutate(date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1")))


# The palette with black:
# cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
#scale_fill_manual(values=cbPalette)
# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

######## Main Island #########
dt <- read.csv("result/oknw_Flu_weather_final/smap_res_3subpref_summary/smap_res_mainisland.csv")
dt$date <- as.Date(dt_ID$date) # date is for time series visualization
dt$year <- dt_ID$year
dt$week <- dt_ID$week
colnames(dt)

date_xaxis <- as.POSIXct(dt$date, format = "%Y-%m-%d")

dt <- dt %>%
  mutate(season = if_else(week < 36, year - 1, year))
dt$season <- as.factor(dt$season)

max_points <- dt %>%
  group_by(season) %>%
  filter(flu_log_ts == max(flu_log_ts, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()

# グラフ作成用の関数を定義
plot_seasonal_timeseries <- function(df, var_name, plot_title) {

  # --- データ加工 ---
  # 関数内で使う変数名を指定するために .data[[var_name]] を使用
  plot_data <- df %>%
    dplyr::group_by(season) %>%
    dplyr::mutate(
      season_avg = mean(.data[[var_name]], na.rm = TRUE)
    ) %>%
    dplyr::ungroup()

  # --- グラフ描画 ---
  p <- ggplot(plot_data, aes(x = date_xaxis, color = factor(season))) +
    # 元データを点線で薄く描画
    geom_line(aes(y = .data[[var_name]])) +
    # シーズンごとの平均値を太い実線で描画
    geom_line(aes(y = season_avg), size = 1) +

    # 共通のグラフ要素
    geom_hline(
      yintercept = quantile(na.omit(df[[var_name]]), c(0.25, 0.5, 0.75)),
      linetype = "dashed",
      color = c("gray", "black", "gray")
    ) +
    geom_vline(
      xintercept = as.POSIXct(c("2012-03-19", "2013-03-11", "2014-03-10", "2015-02-23", "2016-03-21", "2017-03-06", "2018-02-26", "2019-02-25")),
      linetype = "dashed",
      color = "blue"
    ) +
    scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
    scale_colour_viridis_d(name = "Season") +

    # 引数から受け取ったタイトルとy軸ラベルを設定
    labs(y = "Coefficient", title = plot_title) +

    theme(
      text = element_text(size = 20),
      axis.text.x  = element_text(angle = 90),
      axis.title.x = element_blank(),
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      legend.position = "right"
    )
  ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/mainisland_",target_variables,".tiff"),
  units = "in",
  width = 24,
  height = 8,
  dpi = 300,
  compression = 'lzw'
)
  # 作成したggplotオブジェクトを返す
  return(p)
}

# 1. グラフ化したい変数名のリストを作成
target_variables <- c("humidity_cause_flu", "temperature_cause_flu", "flu_log_ts")

# 2. For loopですべての変数に対してグラフを作成し、表示する
for (variable in target_variables) {

  # 関数を呼び出してプロットを作成
  my_plot <- plot_seasonal_timeseries(
    df = dt,
    var_name = variable,
    plot_title = paste("Time series of", variable) # タイトルを自動生成
  )

  # プロットを表示
  print(my_plot)
}

##----##
# 年別（シーズン別）の移動平均と移動中央値を計算
# k=4 は4週間の移動平均/中央値を意味します。必要に応じて変更してください。
dt_moving <- dt %>%
  arrange(season, date) %>% # 各シーズン内で日付順に並べ替え
  group_by(season) %>%
  mutate(
    moving_avg_humid = rollmean(humidity_cause_flu, k = 16, fill = NA, align = "right"),
    moving_median_humid = rollmedian(humidity_cause_flu, k = 16, fill = NA, align = "right"),
    moving_avg_temp = rollmean(temperature_cause_flu, k = 16, fill = NA, align = "right"),
    moving_median_temp = rollmedian(temperature_cause_flu, k = 16, fill = NA, align = "right")
  ) %>%
  ungroup()

# ggplotの描画
humidFlu <- ggplot(dt_moving, aes(x = date_xaxis)) +
  # 元のデータ
  geom_line(aes(y = humidity_cause_flu, color = factor(season)), alpha = 0.5) + # 元の線を少し薄くする

  # 4週移動平均線を追加
  # geom_line(aes(y = moving_avg_humid, color = factor(season)), linetype = "solid", size = 0.5) +

  # 4週移動中央値線を追加 (破線で表示)
  geom_line(aes(y = moving_median_humid, color = factor(season)), linetype = "dashed", size = 0.5) +

  # 元のコードの他の要素
  geom_hline(
    yintercept = quantile(na.omit(dt$humidity_cause_flu), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(
    xintercept = as.POSIXct(c("2012-03-19", "2013-03-11", "2014-03-10", "2015-02-23", "2016-03-21", "2017-03-06", "2018-02-26", "2019-02-25")),
    linetype = "dashed",
    color = "blue"
  ) +
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Coefficient", title = "Humidity cause influenza with Moving Average/Median (4-week)") +
  scale_colour_manual(values = cbPalette) +
  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5),
    legend.position = "none" # 凡例が多いので非表示にする（必要なら "right" などに変更）
  )

# グラフの表示
print(humidFlu)

# 年別（シーズン別）の平均値を計算
dt_seasonal_mean <- dt %>%
  dplyr::group_by(season) %>%
  dplyr::mutate(
    # 各シーズン内のhumidity_cause_fluの平均値を計算
    # na.rm = TRUE は欠損値を無視して計算するオプション
    season_avg_humid = mean(humidity_cause_flu, na.rm = TRUE),
    season_avg_temp = mean(temperature_cause_flu, na.rm = TRUE)
  ) %>%
  ungroup()

# ggplotの描画
humidFlu_mean <- ggplot(dt_seasonal_mean, aes(x = date_xaxis, color = factor(season))) +
  # 元のデータを薄く背景として描画
  geom_line(aes(y = humidity_cause_flu)) +

  # --- ここが変更点 ---
  # シーズンごとの平均値を太い実線で描画
  geom_line(aes(y = season_avg_humid), size = 1.2) +

  # 元のコードの他の要素
  geom_hline(
    yintercept = quantile(na.omit(dt$humidity_cause_flu), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(
    xintercept = as.POSIXct(c("2012-03-19", "2013-03-11", "2014-03-10", "2015-02-23", "2016-03-21", "2017-03-06", "2018-02-26", "2019-02-25")),
    linetype = "dashed",
    color = "blue"
  ) +
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Coefficient", title = "Humidity cause influenza with Seasonal Mean") +

  # 色弱に配慮したviridisカラーパレットを適用
  scale_colour_viridis_d(name = "Season") + # 凡例のタイトルを設定

  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5),
    legend.position = "right" # 凡例を表示
  )

# グラフの表示
print(humidFlu_mean)

tempFlu_mean <- ggplot(dt_seasonal_mean, aes(x = date_xaxis, color = factor(season))) +
  # 元のデータを薄く背景として描画
  geom_line(aes(y = temperature_cause_flu)) +

  # --- ここが変更点 ---
  # シーズンごとの平均値を太い実線で描画
  geom_line(aes(y = season_avg_temp), size = 1.2) +

  # 元のコードの他の要素
  geom_hline(
    yintercept = quantile(na.omit(dt$temperature_cause_flu), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(
    xintercept = as.POSIXct(c("2012-03-19", "2013-03-11", "2014-03-10", "2015-02-23", "2016-03-21", "2017-03-06", "2018-02-26", "2019-02-25")),
    linetype = "dashed",
    color = "blue"
  ) +
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Coefficient", title = "Temperature cause influenza with Seasonal Average, Main Island") +

  # 色弱に配慮したviridisカラーパレットを適用
  scale_colour_viridis_d(name = "Season") + # 凡例のタイトルを設定

  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5),
    legend.position = "right" # 凡例を表示
  )

# グラフの表示
print(humidFlu_mean)

# ggplotの描画
tempFlu <- ggplot(dt_moving, aes(x = date_xaxis)) +
  # 元のデータ
  geom_line(aes(y = temperature_cause_flu, color = factor(season)), alpha = 0.5) + # 元の線を少し薄くする

  # 4週移動平均線を追加
  # geom_line(aes(y = moving_avg_temp, color = factor(season)), linetype = "solid", size = 0.5) +

  # 4週移動中央値線を追加 (破線で表示)
  geom_line(aes(y = moving_median_temp, color = factor(season)), linetype = "dashed", size = 0.5) +

  # 元のコードの他の要素
  geom_hline(
    yintercept = quantile(na.omit(dt$temperature_cause_flu), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(
    xintercept = as.POSIXct(c("2012-03-19", "2013-03-11", "2014-03-10", "2015-02-23", "2016-03-21", "2017-03-06", "2018-02-26", "2019-02-25")),
    linetype = "dashed",
    color = "blue"
  ) +
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Coefficient", title = "Temprature cause influenza with Moving Average/Median (4-week)") +
  scale_colour_manual(values = cbPalette) +
  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5),
    legend.position = "none" # 凡例が多いので非表示にする（必要なら "right" などに変更）
  )

# グラフの表示
print(tempFlu)

##----##

humidFlu <- ggplot(dt, aes(date_xaxis)) +
  geom_line(aes(y = humidity_cause_flu,color = factor(season))) +
  geom_hline(
    yintercept = quantile(na.omit(dt$humidity_cause_flu), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(xintercept = as.POSIXct(c("2012-03-19","2013-03-11","2014-03-10","2015-02-23","2016-03-21","2017-03-06","2018-02-26",
                                       "2019-02-25")),
             linetype = "dashed",
             color = "blue")　+
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Coefficient") +
  ggtitle("Humidity cause influenza") + 
  scale_colour_manual(values=cbPalette)+
  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5)
  )
ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/mainisland_Humidity cause influenza.tiff"),
  units = "in",
  width = 24,
  height = 8,
  dpi = 300,
  compression = 'lzw'
)
####
tempFlu <- ggplot(dt, aes(date_xaxis)) +
  geom_line(aes(y = temperature_cause_flu,color = factor(season))) +
  geom_hline(
    yintercept = quantile(na.omit(dt$temperature_cause_flu), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(xintercept = as.POSIXct(c("2012-03-19","2013-03-11","2014-03-10","2015-02-23","2016-03-21","2017-03-06","2018-02-26",
                                       "2019-02-25")),
             linetype = "dashed",
             color = "blue")　+
  scale_colour_manual(values=cbPalette)+
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Coefficient") +
  ggtitle("Temprature cause influenza") + 
  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5)
  )
ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/mainisland_Temperature cause influenza.tiff"),
  units = "in",
  width = 24,
  height = 8,
  dpi = 300,
  compression = 'lzw'
)
####
logFlu <- ggplot(dt, aes(date_xaxis)) +
  geom_line(aes(y = flu_log_ts,color = factor(season))) +
  geom_hline(
    yintercept = quantile(na.omit(dt$flu_log_ts), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(xintercept = as.POSIXct(c("2012-03-19","2013-03-11","2014-03-10","2015-02-23","2016-03-21","2017-03-06","2018-02-26",
                                       "2019-02-25")),
             linetype = "dashed",
             color = "blue")　+
  scale_colour_manual(values=cbPalette)+
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Incidence") +
  ggtitle("Log-Transformed Influenza Incidence") + 
  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5)
  )
ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/mainisland_Log-Transformed Influenza Incidence.tiff"),
  units = "in",
  width = 24,
  height = 8,
  dpi = 300,
  compression = 'lzw'
)

####
merged_plot <- plot_grid(logFlu,
                         tempFlu,
                         humidFlu,
                         align = "v",
                         ncol = 1)
ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/mainisland_merged.tiff"),
  units = "in",
  width = 36,
  height = 24,
  dpi = 300,
  compression = 'lzw'
)


######## Miyako #########
dt <- read.csv("result/oknw_Flu_weather/smap_res_3subpref_summary/smap_res_miyako.csv")
dt$date <- as.Date(dt_ID$date) # date is for time series visualization
dt$year <- dt_ID$year
dt$week <- dt_ID$week
colnames(dt)

date_xaxis <- as.POSIXct(dt$date, format = "%Y-%m-%d")

dt <- dt %>%
  mutate(season = if_else(week < 36, year - 1, year))

max_points <- dt %>%
  group_by(season) %>%
  filter(flu_log_ts == max(flu_log_ts, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()

##----##

rainFlu <- ggplot(dt, aes(date_xaxis)) +
  geom_line(aes(y = percipitation_cause_flu,color = factor(season))) +
  geom_hline(
    yintercept = quantile(na.omit(dt$percipitation_cause_flu), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(xintercept = as.POSIXct(c("2012-01-30","2013-01-21","2014-01-13","2015-01-12","2016-02-01","2017-02-06","2018-01-15","2019-01-14")),
             linetype = "dashed",
             color = "blue")　+
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Coefficient") +
  ggtitle("Percipitation cause influenza") + 
  scale_colour_manual(values=cbPalette)+
  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5)
  )
ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/miyako_Percipitation cause influenza.tiff"),
  units = "in",
  width = 24,
  height = 8,
  dpi = 300,
  compression = 'lzw'
)
####
tempFlu <- ggplot(dt, aes(date_xaxis)) +
  geom_line(aes(y = temperature_cause_flu,color = factor(season))) +
  geom_hline(
    yintercept = quantile(na.omit(dt$temperature_cause_flu), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(xintercept = as.POSIXct(c("2012-01-30","2013-01-21","2014-01-13","2015-01-12","2016-02-01","2017-02-06","2018-01-15","2019-01-14")),
             linetype = "dashed",
             color = "blue")　+
  scale_colour_manual(values=cbPalette)+
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Coefficient") +
  ggtitle("Temprature cause influenza") + 
  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5)
  )
ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/miyako_Temperature cause influenza.tiff"),
  units = "in",
  width = 24,
  height = 8,
  dpi = 300,
  compression = 'lzw'
)
####
logFlu <- ggplot(dt, aes(date_xaxis)) +
  geom_line(aes(y = flu_log_ts,color = factor(season))) +
  geom_hline(
    yintercept = quantile(na.omit(dt$flu_log_ts), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(xintercept = as.POSIXct(c("2012-01-30","2013-01-21","2014-01-13","2015-01-12","2016-02-01","2017-02-06","2018-01-15","2019-01-14")),
             linetype = "dashed",
             color = "blue")　+
  scale_colour_manual(values=cbPalette)+
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Incidence") +
  ggtitle("Log-Transformed Influenza Incidence") + 
  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5)
  )
ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/miyako_Log-Transformed Influenza Incidence.tiff"),
  units = "in",
  width = 24,
  height = 8,
  dpi = 300,
  compression = 'lzw'
)

####
merged_plot <- plot_grid(logFlu,
                         tempFlu,
                         rainFlu,
                         align = "v",
                         ncol = 1)
ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/miyako_merged.tiff"),
  units = "in",
  width = 36,
  height = 24,
  dpi = 300,
  compression = 'lzw'
)


######## Yaeyama #########
dt <- read.csv("result/oknw_Flu_weather/smap_res_3subpref_summary/smap_res_yaeyama.csv")
dt$date <- as.Date(dt_ID$date) # date is for time series visualization
dt$year <- dt_ID$year
dt$week <- dt_ID$week
colnames(dt)

date_xaxis <- as.POSIXct(dt$date, format = "%Y-%m-%d")

dt <- dt %>%
  mutate(season = if_else(week < 36, year - 1, year))

max_points <- dt %>%
  group_by(season) %>%
  filter(flu_log_ts == max(flu_log_ts, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()

##----##

humidFlu <- ggplot(dt, aes(date_xaxis)) +
  geom_line(aes(y = humidity_cause_flu,color = factor(season))) +
  geom_hline(
    yintercept = quantile(na.omit(dt$humidity_cause_flu), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(xintercept = as.POSIXct(c("2012-02-13","2013-03-04","2014-02-10","2015-02-02","2016-03-07","2017-03-13","2018-01-29","2019-01-21")),
             linetype = "dashed",
             color = "blue")　+
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Coefficient") +
  ggtitle("Humidity cause influenza") + 
  scale_colour_manual(values=cbPalette)+
  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5)
  )
ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/yaeyama_Humidity cause influenza.tiff"),
  units = "in",
  width = 24,
  height = 8,
  dpi = 300,
  compression = 'lzw'
)
####
tempFlu <- ggplot(dt, aes(date_xaxis)) +
  geom_line(aes(y = temperature_cause_flu,color = factor(season))) +
  geom_hline(
    yintercept = quantile(na.omit(dt$temperature_cause_flu), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(xintercept = as.POSIXct(c("2012-02-13","2013-03-04","2014-02-10","2015-02-02","2016-03-07","2017-03-13","2018-01-29","2019-01-21")),
             linetype = "dashed",
             color = "blue")　+
  scale_colour_manual(values=cbPalette)+
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Coefficient") +
  ggtitle("Temprature cause influenza") + 
  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5)
  )
ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/yaeyama_Temperature cause influenza.tiff"),
  units = "in",
  width = 24,
  height = 8,
  dpi = 300,
  compression = 'lzw'
)
####
logFlu <- ggplot(dt, aes(date_xaxis)) +
  geom_line(aes(y = flu_log_ts,color = factor(season))) +
  geom_hline(
    yintercept = quantile(na.omit(dt$flu_log_ts), c(0.25, 0.5, 0.75)),
    linetype = "dashed",
    color = c("gray", "black", "gray")
  ) +
  geom_vline(xintercept = as.POSIXct(c("2012-02-13","2013-03-04","2014-02-10","2015-02-02","2016-03-07","2017-03-13","2018-01-29","2019-01-21")),
             linetype = "dashed",
             color = "blue")　+
  scale_colour_manual(values=cbPalette)+
  scale_x_datetime(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(y = "Incidence") +
  ggtitle("Log-Transformed Influenza Incidence") + 
  theme(
    text = element_text(size = 20),
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5)
  )
ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/yaeyama_Log-Transformed Influenza Incidence.tiff"),
  units = "in",
  width = 24,
  height = 8,
  dpi = 300,
  compression = 'lzw'
)

####
merged_plot <- plot_grid(logFlu,
                         tempFlu,
                         humidFlu,
                         align = "v",
                         ncol = 1)
ggsave(
  paste0("result/oknw_Flu_weather/smap_res_3subpref_summary/figure/yaeyama_merged.tiff"),
  units = "in",
  width = 36,
  height = 24,
  dpi = 300,
  compression = 'lzw'
)