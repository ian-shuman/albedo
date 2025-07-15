###########################################################################################
#
#        This script recreates figure S2 in Shuman et al.
#
#    --- Last updated:  2025.07.14 By Ian Shuman <ins2109@columbia.edu>
###########################################################################################


library(ncdf4) 
library(raster) 
library(rgdal) 
library(ggplot2) 
library(terra)
library(dplyr)
library(stringr)
library(lubridate)
library(zoo)
library(tidyr)

raw_albedo_ts <- rast("~/Downloads/albedo/Data/council_watershed_albedo_ts_orignal.dat")  
smoothed_albedo_ts <- rast("~/Downloads/albedo/Data/council_watershed_albedo_ts_smoothed.dat")  
raw_albedo_ts <- raw_albedo_ts/10000

vals_mat <- as.matrix(raw_albedo_ts)
smoothed_mat <- as.matrix(smoothed_albedo_ts)


# Convert to long format data frame: columns = pixel_id, band/time, value
vals_df <- as.data.frame(vals_mat) %>%
  mutate(pixel_id = row_number()) %>%
  pivot_longer(
    cols = -pixel_id,
    names_to = "time",
    values_to = "value"
  )
smoothed_df <- as.data.frame(smoothed_mat) %>%
  mutate(pixel_id = row_number()) %>%
  pivot_longer(
    cols = -pixel_id,
    names_to = "time",
    values_to = "value"
  )

# Convert 'time' from factor/character like "lyr.1" to numeric band index
vals_df$time <- as.numeric(gsub("lyr.", "", vals_df$time))
vals_df<- vals_df %>% filter(value <= 1.0)
smoothed_df$time <- as.numeric(gsub("lyr.", "", smoothed_df$time))
smoothed_df<- smoothed_df %>% filter(value <= 1.0)

smooth_stats_df <- smoothed_df %>%
  group_by(time) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE)
  )

ggplot() +
  geom_errorbar(
    data = smooth_stats_df,
    aes(x = time, ymin = min, ymax = max, color = "Smoothed + Filtered (Range)"),
    height = 0.1, size = 1, alpha = 0.4,
    inherit.aes = FALSE) +
  geom_point(
    data = vals_df,
    aes(x = time, y = value, color = "Raw Albedo Observations"),
    alpha = 0.1, size = 0.8) +
  geom_point(
    data = smooth_stats_df,
    aes(x = time, y = mean, color = "Smoothed + Filtered (Mean)"),
    size = 1) +
  labs(
    title = "Raw v.s. Smoothed Albedo Time Series",
    x = "DOY", y = "Albedo", color = "Legend") +
  scale_color_manual(values = c(
    "Raw Albedo Observations" = "grey30",
    "Smoothed + Filtered (Mean)" = "blue",
    "Smoothed + Filtered (Range)" = "red")) +
  theme_classic() +
  theme(axis.text = element_text(size=21, color = 'black'),
        axis.title=element_text(size=25),
        axis.text.x = element_text(angle = 60, hjust = 1), 
        aspect.ratio = 1, 
        legend.text = element_text(size = 14))

#Calculate how many points lie outside of the smoothed range
smoothed_min <- global(smoothed_albedo_ts, "min", na.rm = TRUE)[1,1]
smoothed_max <- global(smoothed_albedo_ts, "max", na.rm = TRUE)[1,1]
raw_vals <- values(raw_albedo_ts, mat = TRUE)
raw_vals[raw_vals > 1 | raw_vals < 0] <- NA
smoothed_vals <- values(smoothed_albedo_ts, mat = TRUE)


n_outside <- sum(raw_vals < smoothed_min | raw_vals > smoothed_max, na.rm = TRUE)
n_total <- sum(!is.na(raw_vals))
prop_outside <- n_outside / n_total
cat("Smoothed time series value range:\n")
cat("  Min:", smoothed_min, "\n")
cat("  Max:", smoothed_max, "\n\n")

cat("Raw values outside this range:\n")
cat("  Count:", n_outside, "\n")
cat("  Total valid points:", n_total, "\n")
cat("  Proportion:", round(prop_outside * 100, 2), "%\n")

#Now find out how many pixels would be excluded using reviewer 2's alternative range of 0.15-0.85
n_outside_r2 <- sum(raw_vals < 0.15 | raw_vals > 0.85, na.rm = TRUE)
prop_outside_r2 <- n_outside_r2 / n_total

cat("Raw values outside this range:\n")
cat("  Count:", n_outside_r2, "\n")
cat("  Total valid points:", n_total, "\n")
cat("  Proportion:", round(prop_outside_r2 * 100, 2), "%\n")


#Calculate the mean and standard deviation of the dataset before and after smoothing and filtering
colnames(raw_vals) <- as.numeric(colnames(raw_vals))
colnames(smoothed_vals) <- as.numeric(colnames(smoothed_vals))

raw_data <- raw_vals[, colnames(raw_vals)]
data_smoothed <- smoothed_vals[, colnames(smoothed_vals)]

raw_values <- unlist(raw_data)
raw_values <- raw_values[!is.na(raw_values)]
values_smoothed <- unlist(data_smoothed)
values_smoothed <- values_smoothed[!is.na(values_smoothed)]

raw_mean <- mean(raw_values)
raw_sd   <- sd(raw_values)
mean_smoothed <- mean(values_smoothed)
sd_smoothed   <- sd(values_smoothed)

cat("Raw Mean:", raw_mean, "\n")
cat("Raw SD:", raw_sd, "\n")
cat("Smoothed Mean:", mean_smoothed, "\n")
cat("Smoothed SD:", sd_smoothed, "\n")