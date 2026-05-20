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
  reshape2,
  ISOweek,
  rEDM,
  ggplot2,
  glue,
  stats,
  dplyr,
  gridExtra,
  cowplot,
  rlang,
  writexl,
  openxlsx,
  lattice
)
packageVersion("rEDM")
# library(rEDM)
# vignette("rEDM-tutorial")


## Set path
# setwd handled by edm_code bootstrap
## plot save path
fig_date_path <- "result/figure/oknw_ID/24-09-11/"
hm_path <- paste0(fig_date_path, "subpref_heatmap/")

num <- 1
col_name <- c("naha_dt", "north_dt", "middle_dt", "south_dt", "myk_dt", "yeym_dt", "oknw_dt", "jp_dt")
col_pattern <- c("_Naha", "_North", "_Middle", "_South", "_Miyako", "_Yaeyama", "_Okinawa", "_Japan")
file_name <- c("Xmean", "Xmedian", "Ymean", "Ymedian")
for (i in 1:length(file_name)) {
  for (num in 1:length(col_name)) {
    data <- read.xlsx(paste0("result/analysis_table/subpref/", file_name[i], "/", col_name[num], "_multiSmap_", file_name[i], ".xlsx"))
    rownames(data) <- data$X
    data$X <- NULL
    dt <- data.matrix(data)
    rownames(dt) <- gsub(col_pattern[num], "", rownames(dt))
    colnames(dt) <- gsub(col_pattern[num], "", colnames(dt))

    df <- melt(dt)
    colnames(df) <- c("X", "Y", "Value")
    ggplot(df, aes(x = Y, y = X, fill = Value)) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
      labs(title = col_name[num]) +
      theme_cowplot() +
      theme(text = element_text(size = 20),
            plot.title = element_text(face = "bold")) +
      geom_raster()
    ggsave(paste0(hm_path, "/",file_name[i],"/heatmap", col_pattern[num], ".tiff"), units = "cm",
           width = 25, height = 20, dpi = 300, compression = 'lzw')
  }
}


i <-2
num<-3
data <- read.xlsx(paste0("result/analysis_table/subpref/", file_name[i], "/", col_name[num], "_multiSmap_", file_name[i], ".xlsx"))
rownames(data) <- data$X
data$X <- NULL
dt <- data.matrix(data)
rownames(dt) <- gsub(col_pattern[num], "", rownames(dt))
colnames(dt) <- gsub(col_pattern[num], "", colnames(dt))

df <- melt(dt)
colnames(df) <- c("X", "Y", "Value")
ggplot(df, aes(x = Y, y = X, fill = Value)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(title = col_name[num]) +
  theme_cowplot() +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold")) +
  geom_raster()
ggsave(paste0(hm_path, "heatmap", col_pattern[num], ".tiff"), units = "cm",
       width = 25, height = 20, dpi = 300, compression = 'lzw')

