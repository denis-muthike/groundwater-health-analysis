###Well drawdown and recovery analysis - Sep 2025; Denis Muthike
##Prepared for the manuscript XXXX submitted to AGU Water Resources Research journal

#Load libraries, if not installed, use install.packages("") to install the following packages and load them into your R environment


# Required libraries
# Data manipulation and wrangling
library(data.table)
library(dplyr)
library(lubridate)
library(tidyr)
library(tidyverse)
library(zoo)

# Spatial data and raster handling
library(terra)
library(raster)
library(rasterVis)
library(sf)
library(exactextractr)
library(sp)

# Visualization
library(ggplot2)
library(ggpubr)
library(gridExtra)

# Statistical modeling and analysis
library(extRemes)
library(mgcv)
library(trend)

# Climate indices and drought analysis
library(SCI)
library(SPEI)

# Miscellaneous
library(reshape2)
setwd("/home/denis/projects/drip") #specify path to project data

gwl_all <- read.csv("gwl_all.csv")

# PARAMETERS 
dt_minutes <- 15            # sampling interval (15 min)
min_continuous_minutes <- 15 # minimum contiguous run length to consider (set to 15 min => at least two obs)

# -----------------------------
# 1. Prepare regular 15-min timeline and align data (do NOT interpolate)
# -----------------------------
data_gwl_kalobeyei <- subset(gwl_all, device=="CM.S2B.1.1.0-N0000485")##gwl df from the raw freewave data

Kalobeyei_trend <- data_gwl_kalobeyei %>%
  dplyr::select(timestamp_2, waterlevel) %>% 
  arrange(timestamp_2)


# Show duplicate timestamps in observed data
duplicates <- Kalobeyei_trend %>%
  group_by(timestamp_2) %>%
  dplyr::filter(n() > 1)

if(nrow(duplicates) > 0) {
  cat("\nDuplicate timestamps found in observed data:\n")
  print(head(duplicates))
  cat("Number of duplicate groups:", n_distinct(duplicates$timestamp_2), "\n")
}


# Alternative simpler approach using rolling median
clean_duplicates_simple <- function(observed_data) {
  observed_data %>%
    arrange(timestamp_2) %>%
    group_by(timestamp_2) %>%
    summarize(
      # Strategy 1: Use value closest to rolling median
      waterlevel = {
        if (n() == 1) {
          waterlevel
        } else {
          # Calculate rolling median of surrounding points
          time_idx <- which(observed_data$timestamp_2 == first(timestamp_2))
          window_start <- max(1, time_idx - 6)  # 1.5 hours before (6 * 15min)
          window_end <- min(nrow(observed_data), time_idx + 6)  # 1.5 hours after
          
          surrounding_vals <- observed_data$waterlevel[window_start:window_end]
          surrounding_median <- median(surrounding_vals, na.rm = TRUE)
          
          # Select duplicate value closest to surrounding median
          differences <- abs(waterlevel - surrounding_median)
          waterlevel[which.min(differences)]
        }
      },
      n_duplicates = n(),
      method_used = ifelse(n() > 1, "rolling_median", "single_value")
    ) %>%
    ungroup() %>%
    arrange(timestamp_2)
}

# cleaned data with duplicates removed
Kalobeyei_trend_clean <- clean_duplicates_simple(Kalobeyei_trend)
Kalobeyei_trend_clean <- Kalobeyei_trend_clean %>% dplyr::select(timestamp_2,waterlevel)

# Build full 15-min sequence from min to max
full_time <- tibble(timestamp_2 = seq(
  from = min(Kalobeyei_trend_clean$timestamp_2, na.rm = TRUE),
  to   = max(Kalobeyei_trend_clean$timestamp_2, na.rm = TRUE),
  by   = paste0(dt_minutes, " min")
))

# Left join so missing timestamps become NA
Kalobeyei_full <- full_time %>%
  left_join(Kalobeyei_trend_clean, by = "timestamp_2")

# -----------------------------
# 2. Identify continuous segments OF VALID OBSERVATIONS
#    A continuous segment = consecutive rows with non-NA waterlevel.
#    We will only compute delta_h inside these segments.
# -----------------------------
Kalobeyei_full <- Kalobeyei_full %>%
  mutate(valid = !is.na(waterlevel)) %>%
  # segment starts when valid switches from FALSE -> TRUE
  mutate(segment_start = valid & !lag(valid, default = FALSE)) %>%
  # segment_id increments for each new valid-run; NA for invalid rows
  mutate(segment_id = ifelse(valid, cumsum(segment_start), NA_integer_))

# -----------------------------
# 3. Compute delta_h only WITHIN continuous segments
# -----------------------------
Kalobeyei_events <- Kalobeyei_full %>%
  group_by(segment_id) %>%         # NA segment_id rows will be kept but group will be single-NA group
  mutate(
    wl_prev = lag(waterlevel),
    delta_h = ifelse(!is.na(segment_id), waterlevel - wl_prev, NA_real_)  # delta inside segments only
  ) %>%
  ungroup()

# Remove rows where delta_h is NA (first obs of each segment and any invalid rows)
Kalobeyei_events_valid <- Kalobeyei_events %>%
  dplyr::filter(!is.na(delta_h))

# Add date / month (for aggregation)
Kalobeyei_events_valid <- Kalobeyei_events_valid %>%
  mutate(
    date = as.Date(timestamp_2),
    month = floor_date(timestamp_2, "month")
  )

# -----------------------------
# 4. Compute DAILY drawdown & recovery (mm/day)
#    - recovery: sum of positive deltas
#    - drawdown: sum of absolute negative deltas
# -----------------------------
Kalobeyei_daily <- Kalobeyei_events_valid %>%
  group_by(date) %>%
  summarise(
    recovery_mm = sum(delta_h[delta_h > 0], na.rm = TRUE),
    drawdown_mm = -sum(delta_h[delta_h < 0], na.rm = TRUE),  # positive value
    # bookkeeping
    n_obs = n(),
    .groups = "drop_last"
  ) %>% ungroup() %>% 
  # If no valid observations for that day, the row will be absent; you can add zeros or NA as needed.
  arrange(date)
# -----------------------------
# 5. Compute MONTHLY drawdown & recovery (mm/month)
# -----------------------------
Kalobeyei_monthly <- Kalobeyei_events_valid %>%
  group_by(month) %>%
  summarise(
    recovery_mm = sum(delta_h[delta_h > 0], na.rm = TRUE),
    drawdown_mm = -sum(delta_h[delta_h < 0], na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  arrange(month)

# Convert to meters for plotting (optional)
Kalobeyei_monthly <- Kalobeyei_monthly %>%
  mutate(
    recovery_m = recovery_mm / 1000,
    drawdown_m = drawdown_mm / 1000,
    net_m = recovery_m - drawdown_m
  )

# Save outputs
write.csv(Kalobeyei_daily, "Kalobeyei_daily_drawdown_recovery.csv", row.names = FALSE)
write.csv(Kalobeyei_monthly, "Kalobeyei_monthly_drawdown_recovery.csv", row.names = FALSE)

# -----------------------------
# 6. Plot monthly bars: recovery above zero, drawdown below zero
# -----------------------------
# Prepare plotting dataframe: ensure drawdown is negative for plotting below zero
Kalobeyei_plot_df <- Kalobeyei_monthly %>%
  transmute(
    month = month,
    recovery_m = recovery_m,
    drawdown_m = -drawdown_m  # negative for plotting
  )

# Monthly recovery/recharge bar plot
ggplot(Kalobeyei_plot_df, aes(x = month)) +
  geom_col(aes(y = recovery_m, fill = "Recovery"), alpha = 0.7) +
  geom_col(aes(y = drawdown_m, fill = "Drawdown"), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(
    title = "Monthly Well Drawdown & Recovery",
    subtitle = "Kalobeyei",
    x = "Month",
    y = "Total Water Column Height (m)",
    fill = ""
  ) +
  scale_y_continuous(labels = function(x) abs(x)) + #to keep values on the y axis positive
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

###Loturerei
# -----------------------------
# 1. Prepare regular 15-min timeline and align data (do NOT interpolate)
# -----------------------------
data_gwl_loturerei <- subset(gwl_all, device=="CM.S2B.1.1.0-N0000481")##gwl df from the raw freewave data

loturerei_trend <- data_gwl_loturerei %>%
  dplyr::select(timestamp_2, waterlevel) %>% 
  arrange(timestamp_2)

# Show duplicate timestamps in observed data
duplicates <- loturerei_trend %>%
  group_by(timestamp_2) %>%
  dplyr::filter(n() > 1)

if(nrow(duplicates) > 0) {
  cat("\nDuplicate timestamps found in observed data:\n")
  print(head(duplicates))
  cat("Number of duplicate groups:", n_distinct(duplicates$timestamp_2), "\n")
}

# cleaned data with duplicates removed
loturerei_trend_clean <- clean_duplicates_simple(loturerei_trend)
loturerei_trend_clean <- loturerei_trend_clean %>% dplyr::select(timestamp_2,waterlevel)

# # Compare before and after
cat("=== DUPLICATE CLEANING RESULTS ===\n")
cat("Original data rows:", nrow(loturerei_trend), "\n")
cat("Cleaned data rows:", nrow(loturerei_trend_clean), "\n")
cat("Duplicate groups removed:", nrow(loturerei_trend) - nrow(loturerei_trend_clean), "\n")


# Build full 15-min sequence from min to max
full_time <- tibble(timestamp_2 = seq(
  from = min(loturerei_trend_clean$timestamp_2, na.rm = TRUE),
  to   = max(loturerei_trend_clean$timestamp_2, na.rm = TRUE),
  by   = paste0(dt_minutes, " min")
))

# Left join so missing timestamps become NA
loturerei_full <- full_time %>%
  left_join(loturerei_trend_clean, by = "timestamp_2")

# -----------------------------
# 2. Identify continuous segments OF VALID OBSERVATIONS
#    A continuous segment = consecutive rows with non-NA waterlevel.
#    We will only compute delta_h inside these segments.
# -----------------------------
loturerei_full <- loturerei_full %>%
  mutate(valid = !is.na(waterlevel)) %>%
  # segment starts when valid switches from FALSE -> TRUE
  mutate(segment_start = valid & !lag(valid, default = FALSE)) %>%
  # segment_id increments for each new valid-run; NA for invalid rows
  mutate(segment_id = ifelse(valid, cumsum(segment_start), NA_integer_))

# -----------------------------
# 3. Compute delta_h only WITHIN continuous segments
# -----------------------------
loturerei_events <- loturerei_full %>%
  group_by(segment_id) %>%         # NA segment_id rows will be kept but group will be single-NA group
  mutate(
    wl_prev = lag(waterlevel),
    delta_h = ifelse(!is.na(segment_id), waterlevel - wl_prev, NA_real_)  # delta inside segments only
  ) %>%
  ungroup()

# Remove rows where delta_h is NA (first obs of each segment and any invalid rows)
loturerei_events_valid <- loturerei_events %>%
  dplyr::filter(!is.na(delta_h))

# Add date / month (for aggregation)
loturerei_events_valid <- loturerei_events_valid %>%
  mutate(
    date = as.Date(timestamp_2),
    month = floor_date(timestamp_2, "month")
  )

# -----------------------------
# 4. Compute DAILY drawdown & recovery (mm/day)
#    - recovery: sum of positive deltas
#    - drawdown: sum of absolute negative deltas
# -----------------------------
loturerei_daily <- loturerei_events_valid %>%
  group_by(date) %>%
  summarise(
    recovery_mm = sum(delta_h[delta_h > 0], na.rm = TRUE),
    drawdown_mm = -sum(delta_h[delta_h < 0], na.rm = TRUE),  # positive value
    n_obs = n(),
    .groups = "drop"
  ) %>%
  # If no valid observations for that day, the row will be absent; you can add zeros or NA as needed.
  arrange(date)

# -----------------------------
# 5. Compute MONTHLY drawdown & recovery (mm/month)
# -----------------------------
loturerei_monthly <- loturerei_events_valid %>%
  group_by(month) %>%
  summarise(
    recovery_mm = sum(delta_h[delta_h > 0], na.rm = TRUE),
    drawdown_mm = -sum(delta_h[delta_h < 0], na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  arrange(month)

# Convert to meters for plotting (optional)
loturerei_monthly <- loturerei_monthly %>%
  mutate(
    recovery_m = recovery_mm / 1000,
    drawdown_m = drawdown_mm / 1000,
    net_m = recovery_m - drawdown_m
  )

# -----------------------------
# 6. Plot monthly bars: recovery above zero, drawdown below zero
# -----------------------------
# Prepare plotting dataframe: ensure drawdown is negative for plotting below zero
loturerei_plot_df <- loturerei_monthly %>%
  transmute(
    month = month,
    recovery_m = recovery_m,
    drawdown_m = -drawdown_m  # negative for plotting
  )

loturerei_plot_df <- loturerei_plot_df[loturerei_plot_df$month <= "2024-10-01",]

# Monthly recovery/recharge bar plot
ggplot(loturerei_plot_df, aes(x = month)) +
  geom_col(aes(y = recovery_m, fill = "Recovery"), alpha = 0.7) +
  geom_col(aes(y = drawdown_m, fill = "Drawdown"), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(
    title = "Monthly Well Drawdown & Recovery",
    subtitle = "Loturerei",
    x = "Month",
    y = "Total Water Column Height (m)",
    fill = ""
  ) +
  scale_y_continuous(labels = function(x) abs(x)) + #to keep values on the y axis positive
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )


###Kakemera
# -----------------------------
# 1. Prepare regular 15-min timeline and align data (do NOT interpolate)
# -----------------------------
data_gwl_kakemera <- subset(gwl_all, device=="CM.S2B.1.1.0-N0000491")##gwl df from the raw freewave data

kakemera_trend <- data_gwl_kakemera %>%
  dplyr::select(timestamp_2, waterlevel) %>% 
  arrange(timestamp_2)

kakemera_trend <- na.omit(kakemera_trend) #remove NAs
# kakemera_trend <- kakemera_trend[kakemera_trend$waterlevel >=0,] #remove spurious waterlevel measurements #negative values

# Show duplicate timestamps in observed data
duplicates <- kakemera_trend %>%
  group_by(timestamp_2) %>%
  dplyr::filter(n() > 1)

if(nrow(duplicates) > 0) {
  cat("\nDuplicate timestamps found in observed data:\n")
  print(head(duplicates))
  cat("Number of duplicate groups:", n_distinct(duplicates$timestamp_2), "\n")
}

# cleaned data with duplicates removed
kakemera_trend_clean <- clean_duplicates_simple(kakemera_trend)
kakemera_trend_clean <- kakemera_trend_clean %>% dplyr::select(timestamp_2,waterlevel)

# # Compare before and after
cat("=== DUPLICATE CLEANING RESULTS ===\n")
cat("Original data rows:", nrow(kakemera_trend), "\n")
cat("Cleaned data rows:", nrow(kakemera_trend_clean), "\n")
cat("Duplicate groups removed:", nrow(kakemera_trend) - nrow(kakemera_trend_clean), "\n")


# Build full 15-min sequence from min to max
full_time <- tibble(timestamp_2 = seq(
  from = min(kakemera_trend_clean$timestamp_2, na.rm = TRUE),
  to   = max(kakemera_trend_clean$timestamp_2, na.rm = TRUE),
  by   = paste0(dt_minutes, " min")
))

# Left join so missing timestamps become NA
kakemera_full <- full_time %>%
  left_join(kakemera_trend_clean, by = "timestamp_2")

# -----------------------------
# 2. Identify continuous segments OF VALID OBSERVATIONS
#    A continuous segment = consecutive rows with non-NA waterlevel.
#    We will only compute delta_h inside these segments.
# -----------------------------
kakemera_full <- kakemera_full %>%
  mutate(valid = !is.na(waterlevel)) %>%
  # segment starts when valid switches from FALSE -> TRUE
  mutate(segment_start = valid & !lag(valid, default = FALSE)) %>%
  # segment_id increments for each new valid-run; NA for invalid rows
  mutate(segment_id = ifelse(valid, cumsum(segment_start), NA_integer_))

# -----------------------------
# 3. Compute delta_h only WITHIN continuous segments
# -----------------------------
kakemera_events <- kakemera_full %>%
  group_by(segment_id) %>%         # NA segment_id rows will be kept but group will be single-NA group
  mutate(
    wl_prev = lag(waterlevel),
    delta_h = ifelse(!is.na(segment_id), waterlevel - wl_prev, NA_real_)  # delta inside segments only
  ) %>%
  ungroup()

# Remove rows where delta_h is NA (first obs of each segment and any invalid rows)
kakemera_events_valid <- kakemera_events %>%
  dplyr::filter(!is.na(delta_h))

# Add date / month (for aggregation)
kakemera_events_valid <- kakemera_events_valid %>%
  mutate(
    date = as.Date(timestamp_2),
    month = floor_date(timestamp_2, "month")
  )

# -----------------------------
# 4. Compute DAILY drawdown & recovery (mm/day)
#    - recovery: sum of positive deltas
#    - drawdown: sum of absolute negative deltas
# -----------------------------
kakemera_daily <- kakemera_events_valid %>%
  group_by(date) %>%
  summarise(
    recovery_mm = sum(delta_h[delta_h > 0], na.rm = TRUE),
    drawdown_mm = -sum(delta_h[delta_h < 0], na.rm = TRUE),  # positive value
    n_obs = n(),
    .groups = "drop"
  ) %>%
  # If no valid observations for that day, the row will be absent; you can add zeros or NA as needed.
  arrange(date)

# -----------------------------
# 5. Compute MONTHLY drawdown & recovery (mm/month)
# -----------------------------
kakemera_monthly <- kakemera_events_valid %>%
  group_by(month) %>%
  summarise(
    recovery_mm = sum(delta_h[delta_h > 0], na.rm = TRUE),
    drawdown_mm = -sum(delta_h[delta_h < 0], na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  arrange(month)

# Convert to meters for plotting (optional)
kakemera_monthly <- kakemera_monthly %>%
  mutate(
    recovery_m = recovery_mm / 1000,
    drawdown_m = drawdown_mm / 1000,
    net_m = recovery_m - drawdown_m
  )

# -----------------------------
# 6. Plot monthly bars: recovery above zero, drawdown below zero
# -----------------------------
# Prepare plotting dataframe: ensure drawdown is negative for plotting below zero
kakemera_plot_df <- kakemera_monthly %>%
  transmute(
    month = month,
    recovery_m = recovery_m,
    drawdown_m = -drawdown_m  # negative for plotting
  )

# Monthly recovery/recharge bar plot
ggplot(kakemera_plot_df, aes(x = month)) +
  geom_col(aes(y = recovery_m, fill = "Recovery"), alpha = 0.7) +
  geom_col(aes(y = drawdown_m, fill = "Drawdown"), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(
    title = "Monthly Well Drawdown & Recovery",
    subtitle = "Kakemera",
    x = "Month",
    y = "Total Water Column Height (m)",
    fill = ""
  ) +
  scale_y_continuous(labels = function(x) abs(x)) + #to keep values on the y axis positive
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )


###Nadapal-IRC
# -----------------------------
# 1. Prepare regular 15-min timeline and align data (do NOT interpolate)
# -----------------------------
data_gwl_nadapal <- subset(gwl_all, device=="CM.S2B.1.1.0-N0000486")##gwl df from the raw freewave data

nadapal_trend <- data_gwl_nadapal %>%
  dplyr::select(timestamp_2, waterlevel) %>% 
  arrange(timestamp_2)

nadapal_trend <- na.omit(nadapal_trend) #remove NAs
# nadapal_trend <- nadapal_trend[nadapal_trend$waterlevel >=0,] #remove spurious waterlevel measurements #negative values

# Show duplicate timestamps in observed data
duplicates <- nadapal_trend %>%
  group_by(timestamp_2) %>%
  dplyr::filter(n() > 1)

if(nrow(duplicates) > 0) {
  cat("\nDuplicate timestamps found in observed data:\n")
  print(head(duplicates))
  cat("Number of duplicate groups:", n_distinct(duplicates$timestamp_2), "\n")
}

# cleaned data with duplicates removed
nadapal_trend_clean <- clean_duplicates_simple(nadapal_trend)
nadapal_trend_clean <- nadapal_trend_clean %>% dplyr::select(timestamp_2,waterlevel)


# Build full 15-min sequence from min to max
full_time <- tibble(timestamp_2 = seq(
  from = floor_date(min(nadapal_trend_clean$timestamp_2, na.rm = TRUE), "15 minutes"),
  to   = ceiling_date(max(nadapal_trend_clean$timestamp_2, na.rm = TRUE), "15 minutes"),
  by   = paste0(dt_minutes, " min")
))

# Left join so missing timestamps become NA
nadapal_full <- full_time %>%
  left_join(nadapal_trend_clean, by = "timestamp_2")

# -----------------------------
# 2. Identify continuous segments OF VALID OBSERVATIONS
#    A continuous segment = consecutive rows with non-NA waterlevel.
#    We will only compute delta_h inside these segments.
# -----------------------------
nadapal_full <- nadapal_full %>%
  mutate(valid = !is.na(waterlevel)) %>%
  # segment starts when valid switches from FALSE -> TRUE
  mutate(segment_start = valid & !lag(valid, default = FALSE)) %>%
  # segment_id increments for each new valid-run; NA for invalid rows
  mutate(segment_id = ifelse(valid, cumsum(segment_start), NA_integer_))

# -----------------------------
# 3. Compute delta_h only WITHIN continuous segments
# -----------------------------
nadapal_events <- nadapal_full %>%
  group_by(segment_id) %>%         # NA segment_id rows will be kept but group will be single-NA group
  mutate(
    wl_prev = lag(waterlevel),
    delta_h = ifelse(!is.na(segment_id), waterlevel - wl_prev, NA_real_)  # delta inside segments only
  ) %>%
  ungroup()

# Remove rows where delta_h is NA (first obs of each segment and any invalid rows)
nadapal_events_valid <- nadapal_events %>%
  dplyr::filter(!is.na(delta_h))

# Add date / month (for aggregation)
nadapal_events_valid <- nadapal_events_valid %>%
  mutate(
    date = as.Date(timestamp_2),
    month = floor_date(timestamp_2, "month")
  )

# -----------------------------
# 4. Compute DAILY drawdown & recovery (mm/day)
#    - recovery: sum of positive deltas
#    - drawdown: sum of absolute negative deltas
# -----------------------------
nadapal_daily <- nadapal_events_valid %>%
  group_by(date) %>%
  summarise(
    recovery_mm = sum(delta_h[delta_h > 0], na.rm = TRUE),
    drawdown_mm = -sum(delta_h[delta_h < 0], na.rm = TRUE),  # positive value
    n_obs = n(),
    .groups = "drop"
  ) %>%
  # If no valid observations for that day, the row will be absent; you can add zeros or NA as needed.
  arrange(date)

# -----------------------------
# 5. Compute MONTHLY drawdown & recovery (mm/month)
# -----------------------------
nadapal_monthly <- nadapal_events_valid %>%
  group_by(month) %>%
  summarise(
    recovery_mm = sum(delta_h[delta_h > 0], na.rm = TRUE),
    drawdown_mm = -sum(delta_h[delta_h < 0], na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  arrange(month)

# Convert to meters for plotting (optional)
nadapal_monthly <- nadapal_monthly %>%
  mutate(
    recovery_m = recovery_mm / 1000,
    drawdown_m = drawdown_mm / 1000,
    net_m = recovery_m - drawdown_m
  )

# -----------------------------
# 6. Plot monthly bars: recovery above zero, drawdown below zero
# -----------------------------
# Prepare plotting dataframe: ensure drawdown is negative for plotting below zero
nadapal_plot_df <- nadapal_monthly %>%
  transmute(
    month = month,
    recovery_m = recovery_m,
    drawdown_m = -drawdown_m  # negative for plotting
  )

# Monthly recovery/recharge bar plot
ggplot(nadapal_plot_df, aes(x = month)) +
  geom_col(aes(y = recovery_m, fill = "Recovery"), alpha = 0.7) +
  geom_col(aes(y = drawdown_m, fill = "Drawdown"), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(
    title = "Monthly Well Drawdown & Recovery",
    subtitle = "Nadapal-IRC",
    x = "Month",
    y = "Total Water Column Height (m)",
    fill = ""
  ) +
  scale_y_continuous(labels = function(x) abs(x)) + #to keep values on the y axis positive
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# End

####COMBINED PLOTS
nadapal_plot_df$site="Nadapal-IRC"
Kalobeyei_plot_df$site = "Kalobeyei"
loturerei_plot_df$site = "Loturerei"
kakemera_plot_df$site = "Kakemera"


all_mon_estimates <- rbind(Kalobeyei_plot_df,loturerei_plot_df,
                           kakemera_plot_df,nadapal_plot_df)

##plot
all_mon_recovery_estimates = ggplot(all_mon_estimates, aes(x = month)) +
  geom_col(aes(y = recovery_m, fill = "Recovery"), alpha = 0.7) +
  geom_col(aes(y = drawdown_m, fill = "Drawdown"), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(x = "",
       y = "Total Water Column Height (m)",
       fill = ""
  ) +facet_wrap(~site, scales="free_y")+#ylim(0, 600)+
  scale_y_continuous(labels = function(x) abs(x)) + #to keep values on the y axis positive
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    axis.title.y = element_text(size = 10)
  )

ggsave("./plots/final-paper/drawdown-recovery-trends-final.pdf", all_mon_recovery_estimates,
       device = cairo_pdf,      # Use ragg device
       width = 7, height = 5, units = "in")           # No scaling


###all sites raw TS
five_days_data_loturerei <- five_days_data_loturerei %>% mutate(site = "Loturerei")
five_days_data_kalobeyei <- five_days_data_kalobeyei %>% mutate(site = "Kalobeyei")
five_days_data_kakemera <- five_days_data_kakemera %>% mutate(site = "Kakemera")
five_days_data_nadapal <- five_days_data_nadapal  %>% mutate(site = "Nadapal-IRC")

all_five_day_ts <- bind_rows(five_days_data_loturerei,five_days_data_kalobeyei,five_days_data_kakemera,five_days_data_nadapal)

y_rng <- range(all_five_day_ts$waterlevel, na.rm = TRUE)

all_five_day_ts_plot = ggplot(all_five_day_ts, aes(x=timestamp_2, y=waterlevel/1000))+
  geom_line()+scale_x_datetime(
    date_labels = "%m/%d %H:%M",  # Includes seconds
    date_breaks = "12 hours"  # Adjust break frequency
  ) +xlab("Time (hourly)")+ylab("Water Column Height (m)")+
  facet_wrap(~ site, scales = "free_x", ncol = 2)+
  scale_y_continuous(limits = y_rng/1000, expand = expansion(mult = c(0, 0.03))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("./plots/final-paper/f_v2/four-site-waterlevel_ts-final.png", plot = all_five_day_ts_plot, device = "png", dpi = 600,
       width = 5, height = 6, units = "in")



#######################################
########Net-storage-analysis
#######################################

#Nadapal
ex_df_daily <- nadapal_events_valid %>%
  group_by(date) %>%
  summarise(
    recovery_mm = sum(delta_h[delta_h > 0], na.rm = TRUE),
    drawdown_mm = -sum(delta_h[delta_h < 0], na.rm = TRUE),  # positive value
    n_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(is_valid = n_obs >= 0.9 * 96) %>%  # check completeness
  filter(is_valid) %>%                      # keep only valid days
  arrange(date)

nadapal_net_storage <- ex_df_daily

# calculate number of days in the month
ex_df_monthly <- ex_df_daily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(
    recovery_mm = sum(recovery_mm, na.rm = TRUE),
    drawdown_mm = sum(drawdown_mm, na.rm = TRUE),
    n_days = n(),  # valid days in this month
    .groups = "drop"
  ) %>%
  mutate(
    days_in_month = days_in_month(month),
    valid_month = n_days >= 0.7 * days_in_month
  ) %>%
  filter(valid_month) %>%
  arrange(month) %>%
  mutate(
    recovery_m = recovery_mm / 1000,
    drawdown_m = drawdown_mm / 1000,
    net_m = recovery_m - drawdown_m
  )

# -----------------------------
# Prepare plotting dataframe: ensure drawdown is negative for plotting below zero
ex_df_plot_df <- ex_df_monthly %>%
  transmute(
    month = month,
    recovery_m = recovery_m,
    drawdown_m = -drawdown_m  # negative for plotting
  )

ex_df_plot_daily <- ex_df_daily %>%
  transmute(
    date = date,
    recovery_m = recovery_mm/1000,
    drawdown_m = -drawdown_mm/1000  # negative for plotting
  )

# Monthly recovery/recharge bar plot
nadapal_daily_estim_plot = ggplot(ex_df_plot_daily, aes(x = date)) +
  geom_col(aes(y = recovery_m, fill = "Recovery"), alpha = 0.7) +
  geom_col(aes(y = drawdown_m, fill = "Drawdown"), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(
    title = "Monthly Well Drawdown & Recovery",
    subtitle = "ex_df-IRC",
    x = "Month",
    y = "Total Water Column Height (m)",
    fill = ""
  ) +
  scale_y_continuous(labels = function(x) abs(x)) + #to keep values on the y axis positive
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

nadapal_daily_estimate <- ex_df_plot_daily

#kalobeyei

ex_df_daily <- Kalobeyei_events_valid %>%
  group_by(date) %>%
  summarise(
    recovery_mm = sum(delta_h[delta_h > 0], na.rm = TRUE),
    drawdown_mm = -sum(delta_h[delta_h < 0], na.rm = TRUE),  # positive value
    n_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(is_valid = n_obs >= 0.9 * 96) %>%  # check completeness
  filter(is_valid) %>%                      # keep only valid days
  arrange(date)


ex_df_monthly <- ex_df_daily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(
    recovery_mm = sum(recovery_mm, na.rm = TRUE),
    drawdown_mm = sum(drawdown_mm, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  ) %>%
  arrange(month)

# Convert to meters for plotting (optional)
ex_df_monthly <- ex_df_monthly %>%
  mutate(
    recovery_m = recovery_mm / 1000,
    drawdown_m = drawdown_mm / 1000,
    net_m = recovery_m - drawdown_m
  )

ex_df_plot_daily <- ex_df_daily %>%
  transmute(
    date = date,
    recovery_m = recovery_mm/1000,
    drawdown_m = -drawdown_mm/1000  # negative for plotting
  )

# -----------------------------
# Prepare plotting dataframe: ensure drawdown is negative for plotting below zero
ex_df_plot_df <- ex_df_monthly %>%
  transmute(
    month = month,
    recovery_m = recovery_m,
    drawdown_m = -drawdown_m  # negative for plotting
  )

# Monthly recovery/recharge bar plot
ggplot(ex_df_plot_daily, aes(x = date)) +
  geom_col(aes(y = recovery_m, fill = "Recovery"), alpha = 0.7) +
  geom_col(aes(y = drawdown_m, fill = "Drawdown"), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(
    title = "Monthly Well Drawdown & Recovery",
    subtitle = "ex_df-IRC",
    x = "Month",
    y = "Total Water Column Height (m)",
    fill = ""
  ) +
  scale_y_continuous(labels = function(x) abs(x)) + #to keep values on the y axis positive
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

kalobeyei_daily_estimate <- ex_df_plot_daily

##Loturerei
ex_df_daily <- loturerei_events_valid %>%
  group_by(date) %>%
  summarise(
    recovery_mm = sum(delta_h[delta_h > 0], na.rm = TRUE),
    drawdown_mm = -sum(delta_h[delta_h < 0], na.rm = TRUE),  # positive value
    n_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(is_valid = n_obs >= 0.9 * 96) %>%  # check completeness
  filter(is_valid) %>%                      # keep only valid days
  arrange(date)

loturerei_net_storage <- ex_df_daily

ex_df_monthly <- ex_df_daily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(
    recovery_mm = sum(recovery_mm, na.rm = TRUE),
    drawdown_mm = sum(drawdown_mm, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  ) %>%
  arrange(month)

# Convert to meters for plotting (optional)
ex_df_monthly <- ex_df_monthly %>%
  mutate(
    recovery_m = recovery_mm / 1000,
    drawdown_m = drawdown_mm / 1000,
    net_m = recovery_m - drawdown_m
  )

ex_df_plot_daily <- ex_df_daily %>%
  transmute(
    date = date,
    recovery_m = recovery_mm/1000,
    drawdown_m = -drawdown_mm/1000  # negative for plotting
  )

# -----------------------------
# Prepare plotting dataframe: ensure drawdown is negative for plotting below zero
ex_df_plot_df <- ex_df_monthly %>%
  transmute(
    month = month,
    recovery_m = recovery_m,
    drawdown_m = -drawdown_m  # negative for plotting
  )

# Monthly recovery/recharge bar plot
ggplot(ex_df_plot_daily, aes(x = date)) +
  geom_col(aes(y = recovery_m, fill = "Recovery"), alpha = 0.7) +
  geom_col(aes(y = drawdown_m, fill = "Drawdown"), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(
    title = "Monthly Well Drawdown & Recovery",
    subtitle = "ex_df-IRC",
    x = "Month",
    y = "Total Water Column Height (m)",
    fill = ""
  ) +
  scale_y_continuous(labels = function(x) abs(x)) + #to keep values on the y axis positive
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

loturerei_daily_estimate <- ex_df_plot_daily

##Kakemera
ex_df_daily <- kakemera_events_valid %>%
  group_by(date) %>%
  summarise(
    recovery_mm = sum(delta_h[delta_h > 0], na.rm = TRUE),
    drawdown_mm = -sum(delta_h[delta_h < 0], na.rm = TRUE),  # positive value
    n_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(is_valid = n_obs >= 0.9 * 96) %>%  # check completeness
  filter(is_valid) %>%                      # keep only valid days
  arrange(date)

kakemera_net_storage <- ex_df_daily

ex_df_monthly <- ex_df_daily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(
    recovery_mm = sum(recovery_mm, na.rm = TRUE),
    drawdown_mm = sum(drawdown_mm, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  ) %>%
  arrange(month)

# Convert to meters for plotting 
ex_df_monthly <- ex_df_monthly %>%
  mutate(
    recovery_m = recovery_mm / 1000,
    drawdown_m = drawdown_mm / 1000,
    net_m = recovery_m - drawdown_m
  )

# -----------------------------
# Prepare plotting dataframe: ensure drawdown is negative for plotting below zero
ex_df_plot_df <- ex_df_monthly %>%
  transmute(
    month = month,
    recovery_m = recovery_m,
    drawdown_m = -drawdown_m  # negative for plotting
  )

ex_df_plot_daily <- ex_df_daily %>%
  transmute(
    date = date,
    recovery_m = recovery_mm/1000,
    drawdown_m = -drawdown_mm/1000  # negative for plotting
  )


# Monthly recovery/recharge bar plot
kakemera_daily_estim_plot = ggplot(ex_df_plot_daily, aes(x = date)) +
  geom_col(aes(y = recovery_m, fill = "Recovery"), alpha = 0.7) +
  geom_col(aes(y = drawdown_m, fill = "Drawdown"), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(
    title = "Monthly Well Drawdown & Recovery",
    subtitle = "ex_df-IRC",
    x = "Month",
    y = "Total Water Column Height (m)",
    fill = ""
  ) +
  scale_y_continuous(labels = function(x) abs(x)) + #to keep values on the y axis positive
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

kakemera_daily_estimate <-ex_df_plot_daily 

###Combined daily-estimates
kalobeyei_daily_estimate <- kalobeyei_daily_estimate %>% mutate(site = "Kalobeyei")
kakemera_daily_estimate <- kakemera_daily_estimate %>% mutate(site = "Kakemera")
nadapal_daily_estimate <- nadapal_daily_estimate %>% mutate(site = "Nadapal-IRC")
loturerei_daily_estimate <- loturerei_daily_estimate %>% mutate(site = "Loturerei")

all_sites <- bind_rows(kalobeyei_daily_estimate, kakemera_daily_estimate, nadapal_daily_estimate, loturerei_daily_estimate)

all_sites_drawdown_recovery = ggplot(all_sites, aes(x = date)) +
  geom_col(aes(y = recovery_m, fill = "Recovery"), alpha = 0.7) +
  geom_col(aes(y = drawdown_m, fill = "Drawdown"), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(x="",y = "Daily total water \n column height (m)",
       fill = ""
  ) +
  scale_y_continuous(labels = function(x) abs(x)) + #to keep values on the y axis positive
  facet_wrap(~ site, scales = "free_y") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# Summarize to site-level means
site_summary_daily <- all_sites %>%
  group_by(site) %>%
  summarise(
    mean_recovery = mean(recovery_m, na.rm = TRUE),
    mean_drawdown = mean((drawdown_m*-1), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Pivot into long format for ggplot
  pivot_longer(cols = starts_with("mean_"),
               names_to = "type",
               values_to = "mean_value") %>%
  mutate(
    type = recode(type,
                  "mean_recovery" = "Recovery",
                  "mean_drawdown" = "Drawdown")
  )

site_summary_month <- all_sites %>%
  mutate(month = month(date)) %>% 
  group_by(site,month) %>%
  summarise(
    mean_recovery = mean(recovery_m, na.rm = TRUE),
    mean_drawdown = mean((drawdown_m*-1), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Pivot into long format for ggplot
  pivot_longer(cols = starts_with("mean_"),
               names_to = "type",
               values_to = "mean_value") %>%
  mutate(
    type = recode(type,
                  "mean_recovery" = "Recovery",
                  "mean_drawdown" = "Drawdown")
  )

all_sites_net_storage = ggplot(site_summary_daily, aes(x = site, y = mean_value, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = "", y = "Mean daily water \n column height (m)",
       fill = ""
  ) +
  theme_minimal() +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12)
  )

ggsave("./plots/final-paper/all_sites_net_storage.png", plot = all_sites_net_storage,
       device = png,
       dpi = 300,
       width = 6, height = 3, units = "in")           # No scaling

ggsave("./plots/final-paper/f_v2/gwl_all_sites_drawdown_recovery.png", plot = all_sites_drawdown_recovery,
       device = png,
       dpi = 600,
       width = 6, height = 4, units = "in")           # No scaling


###weekly anomalies - drawdown and recovery
# Reshape into long format: one column for values, one for category
weekly_df_long <- all_sites
weekly_df_long <- weekly_df_long %>% mutate(drawdown_m = drawdown_m*-1)

weekly_df_long <- weekly_df_long %>%
  pivot_longer(cols = c(recovery_m, drawdown_m),
               names_to = "category",
               values_to = "value")

# Step 1: Aggregate to weekly means (per site and category)
df_weekly <- weekly_df_long %>%
  mutate(week = isoweek(date),
         year = year(date)) %>%
  group_by(site, category, year, week) %>%
  summarise(weekly_mean = mean(value, na.rm = TRUE),
            weekly_total = sum(value, na.rm = TRUE),
            .groups = "drop")

# Step 2: Compute long-term weekly climatology for each site + category
df_weekly <- df_weekly %>%
  group_by(site, category, week) %>%
  mutate(climatology = mean(weekly_mean, na.rm = TRUE)) %>%
  ungroup()

# Step 3: Compute anomalies
df_weekly <- df_weekly %>%
  mutate(anomaly = weekly_mean - climatology)

#plot
##non-stacked bars
ggplot(df_weekly, aes(
  x = as.Date(paste(year, week, 1, sep = "-"), "%Y-%U-%u"),
  y = anomaly,
  fill = category
))+ geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_col(position = position_dodge(width = 12),width = 12) +   # side-by-side bars
  facet_wrap(~site, scales = "free_y") +
  labs(
    x = "",
    y = "Weekly anomaly (m)",
    title = "Weekly recovery and drawdown anomalies by site"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("drawdown_m" = "red", "recovery_m" = "blue"))

##just recovery
df_weekly_recovery <- df_weekly %>% subset(category =="recovery_m")

ggplot(df_weekly_recovery, aes(
  x = as.Date(paste(year, week, 1, sep = "-"), "%Y-%U-%u"),
  y = anomaly,
  fill = category
))+ geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_col(position = position_dodge(width = 12),width = 12) +   # side-by-side bars
  facet_wrap(~site, scales = "free_y") +
  labs(
    x = "",
    y = "Weekly anomaly (m)",
    title = ""
  ) +
  theme_bw() +
  scale_fill_manual(values = c("recovery_m" = "blue"))


###weekly mean d/r
##option 

df_weekly_plot <- df_weekly %>%
  mutate(adj_mean_value = ifelse(category == "drawdown_m", 
                                 -abs(weekly_mean),  # Force negative for drawdown
                                 abs(weekly_mean))) %>%  # Force positive for recovery
  mutate(adj_tot_value = ifelse(category == "drawdown_m", 
                                -abs(weekly_total),  # Force negative for drawdown
                                abs(weekly_total))) %>%  # Force positive for recovery
  mutate(category_z = ifelse(category == "drawdown_m","Drawdown","Recovery"))

##plot
weekly_dd_rvy = ggplot(df_weekly_plot, aes(x = as.Date(paste(year, week, 1, sep="-"), "%Y-%U-%u"),
                                           y = adj_mean_value,
                                           fill = category_z)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_col(position = "identity",alpha = 1) +
  facet_wrap(~site, scales = "free_y") +
  labs(x = "",y = "Weekly mean water column height (m)",
       title = "") +
  scale_y_continuous(labels = function(x) abs(x)) + #to keep values on the y axis positive
  theme_bw() +
  scale_fill_manual(values = c("Drawdown" = "brown", "Recovery" = "blue"))+
  theme(legend.title = element_blank(),
        legend.position = "bottom")  # Add some padding

ggsave("./plots/final-paper/f_v2/weekly-drawdown-recovery-final.png", plot = weekly_dd_rvy, device = "png", dpi = 300,
       width = 6, height = 4, units = "in")


###number of observation months for each site

df_months <- all_sites %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(site) %>%
  summarise(
    n_months = n_distinct(month),
    n_days = n_distinct(date),
    start_date = min(date),
    end_date = max(date)
  )
print(df_months)


####################################
##Groundwater head-rainfall analysis
##

kalobeyei_df <- Kalobeyei_full ##carried over from the drawdown & recovery analysis

# Convert timestamp
kalobeyei_df <- kalobeyei_df %>% dplyr::select(timestamp_2,waterlevel)

# Sort by time
kalobeyei_df <- kalobeyei_df %>% arrange(timestamp_2)

# 1. Identify gaps - confirm gaps filled?
kalobeyei_df <- kalobeyei_df %>%
  mutate(dt_min = as.numeric(difftime(timestamp_2, lag(timestamp_2), units = "mins")),
         gap_flag = ifelse(!is.na(dt_min) & dt_min > 30, TRUE, FALSE))  # >30 min gap

# 2. Compute rate of change
kalobeyei_df <- kalobeyei_df %>%
  mutate(dlevel = waterlevel - lag(waterlevel),
         dt_hr = as.numeric(difftime(timestamp_2, lag(timestamp_2), units = "hours")),
         rate = dlevel / dt_hr)

# 3. Define recovery (non-pumping) = positive slope or near-zero slope
kalobeyei_df <- kalobeyei_df %>%
  mutate(is_recovery = ifelse(!is.na(rate) & rate >= 0, TRUE, FALSE))

# 4. Keep recovery-only data
kalobeyei_df_recovery <- kalobeyei_df %>% dplyr::filter(is_recovery)

# 5. Compute daily 95th percentile from recovery data only
kalobeyei_daily_95th <- kalobeyei_df_recovery %>%
  group_by(date = as.Date(timestamp_2)) %>%
  summarise(
    p95_head = quantile(waterlevel, 0.95, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_obs > 5 &  p95_head  > 100)  # require at least 5 recovery obs per day and remove low values that seem suspect

# 6. Plot
ggplot(kalobeyei_daily_95th, aes(x=date, y=p95_head/1000)) + #convert units to m
  geom_line(color="darkblue") +
  geom_point(color="steelblue", size=0.8) +
  labs(title="Kalobeyei",
       x="Date", y="95th percentile water column height (m)") +
  theme_minimal()


##monthly
kalobeyei_monthly_95th <- kalobeyei_daily_95th %>% 
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>% 
  summarise(mon_head = mean(p95_head, na.rm=T))

ggplot(kalobeyei_monthly_95th, aes(x=month, y=mon_head/1000)) +
  geom_line(color="darkblue") +
  geom_point(color="steelblue", size=0.8) +
  labs(title="Kalobeyei",
       x="Month", y="Mean of daily 95th percentile water column height (m)") +
  theme_minimal()

##Loturerei
loturerei_df <- loturerei_full

# Convert timestamp
loturerei_df <- loturerei_df %>% dplyr::select(timestamp_2,waterlevel)

# Sort by time
loturerei_df <- loturerei_df %>% arrange(timestamp_2)

# 1. Identify gaps
loturerei_df <- loturerei_df %>%
  mutate(dt_min = as.numeric(difftime(timestamp_2, lag(timestamp_2), units = "mins")),
         gap_flag = ifelse(!is.na(dt_min) & dt_min > 30, TRUE, FALSE))  # >30 min gap

# 2. Compute rate of change
loturerei_df <- loturerei_df %>%
  mutate(dlevel = waterlevel - lag(waterlevel),
         dt_hr = as.numeric(difftime(timestamp_2, lag(timestamp_2), units = "hours")),
         rate = dlevel / dt_hr)

# 3. Define recovery (non-pumping) = positive slope or near-zero slope
loturerei_df <- loturerei_df %>%
  mutate(is_recovery = ifelse(!is.na(rate) & rate >= 0, TRUE, FALSE))

# 4. Keep recovery-only data
loturerei_df_recovery <- loturerei_df %>% dplyr::filter(is_recovery)

# 5. Compute daily 95th percentile from recovery data only
loturerei_daily_95th <- loturerei_df_recovery %>%
  group_by(date = as.Date(timestamp_2)) %>%
  summarise(
    p95_head = quantile(waterlevel, 0.95, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_obs > 5 &  p95_head  > 100)  # require at least 5 recovery obs per day

# 6. Plot
ggplot(loturerei_daily_95th, aes(x=date, y=p95_head/1000)) +
  geom_line(color="darkblue") +
  geom_point(color="steelblue", size=0.8) +
  labs(title="Loturerei",
       x="Date", y="95th percentile water column height (m)") +
  theme_minimal()

##monthly
loturerei_monthly_95th <- loturerei_daily_95th %>% 
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>% 
  summarise(mon_head = mean(p95_head, na.rm=T))

ggplot(loturerei_monthly_95th, aes(x=month, y=mon_head/1000)) +
  geom_line(color="darkblue") +
  geom_point(color="steelblue", size=0.8) +
  labs(title="Loturerei",
       x="Month", y="Mean of daily 95th percentile water column height (m)") +
  theme_minimal()


##Kakemera
kakemera_df <- kakemera_full

# Convert timestamp
kakemera_df <- kakemera_df %>% dplyr::select(timestamp_2,waterlevel)

# Sort by time
kakemera_df <- kakemera_df %>% arrange(timestamp_2)

# 1. Identify gaps
kakemera_df <- kakemera_df %>%
  mutate(dt_min = as.numeric(difftime(timestamp_2, lag(timestamp_2), units = "mins")),
         gap_flag = ifelse(!is.na(dt_min) & dt_min > 30, TRUE, FALSE))  # >30 min gap

# 2. Compute rate of change
kakemera_df <- kakemera_df %>%
  mutate(dlevel = waterlevel - lag(waterlevel),
         dt_hr = as.numeric(difftime(timestamp_2, lag(timestamp_2), units = "hours")),
         rate = dlevel / dt_hr)

# 3. Define recovery (non-pumping) = positive slope or near-zero slope
kakemera_df <- kakemera_df %>%
  mutate(is_recovery = ifelse(!is.na(rate) & rate >= 0, TRUE, FALSE))

# 4. Keep recovery-only data
kakemera_df_recovery <- kakemera_df %>% dplyr::filter(is_recovery)

# 5. Compute daily 95th percentile from recovery data only
kakemera_daily_95th <- kakemera_df_recovery %>%
  group_by(date = as.Date(timestamp_2)) %>%
  summarise(
    p95_head = quantile(waterlevel, 0.95, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_obs > 5 &  p95_head  > 100)  # require at least 5 recovery obs per day

# 6. Plot
ggplot(kakemera_daily_95th, aes(x=date, y=p95_head/1000)) +
  geom_line(color="darkblue") +
  geom_point(color="steelblue", size=0.8) +
  labs(title="Kakemera",
       x="Date", y="95th percentile water column height (m)") +
  theme_minimal()

##monthly
kakemera_monthly_95th <- kakemera_daily_95th %>% 
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>% 
  summarise(mon_head = mean(p95_head, na.rm=T))

ggplot(kakemera_monthly_95th, aes(x=month, y=mon_head/1000)) +
  geom_line(color="darkblue") +
  geom_point(color="steelblue", size=0.8) +
  labs(title="Kakemera",
       x="Month", y="Mean of daily 95th percentile water column height (m)") +
  theme_minimal()



##Nadapal
nadapal_df <- nadapal_full

# Convert timestamp
nadapal_df <- nadapal_df %>% dplyr::select(timestamp_2,waterlevel)

# Sort by time
nadapal_df <- nadapal_df %>% arrange(timestamp_2)

# 1. Identify gaps
nadapal_df <- nadapal_df %>%
  mutate(dt_min = as.numeric(difftime(timestamp_2, lag(timestamp_2), units = "mins")),
         gap_flag = ifelse(!is.na(dt_min) & dt_min > 30, TRUE, FALSE))  # >30 min gap

# 2. Compute rate of change
nadapal_df <- nadapal_df %>%
  mutate(dlevel = waterlevel - lag(waterlevel),
         dt_hr = as.numeric(difftime(timestamp_2, lag(timestamp_2), units = "hours")),
         rate = dlevel / dt_hr)

# 3. Define recovery (non-pumping) = positive slope or near-zero slope
nadapal_df <- nadapal_df %>%
  mutate(is_recovery = ifelse(!is.na(rate) & rate >= 0, TRUE, FALSE))

# 4. Keep recovery-only data
nadapal_df_recovery <- nadapal_df %>% dplyr::filter(is_recovery)

# 5. Compute daily 95th percentile from recovery data only
nadapal_daily_95th <- nadapal_df_recovery %>%
  group_by(date = as.Date(timestamp_2)) %>%
  summarise(
    p95_head = quantile(waterlevel, 0.95, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_obs > 5 &  p95_head  > 100)  # require at least 5 recovery obs per day

# 6. Plot
ggplot(nadapal_daily_95th, aes(x=date, y=p95_head/1000)) +
  geom_line(color="darkblue") +
  geom_point(color="steelblue", size=0.8) +
  labs(title="Borehole-9",
       x="Date", y="95th percentile of daily water column height (m)") +
  theme_minimal()

nadapal_monthly_95th <- nadapal_daily_95th %>% 
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>% 
  summarise(mon_head = mean(p95_head, na.rm=T))

ggplot(nadapal_monthly_95th, aes(x=month, y=mon_head/1000)) +
  geom_line(color="darkblue") +
  geom_point(color="steelblue", size=0.8) +
  labs(title="Borehole-9",
       x="Month", y="Mean of daily 95th percentile water column height (m)") +
  theme_minimal()
##

nadapal_daily_95th <- nadapal_daily_95th %>% mutate(site = "Nadapal-IRC")
kakemera_daily_95th <- kakemera_daily_95th %>% mutate(site = "Kakemera")
kalobeyei_daily_95th <- kalobeyei_daily_95th %>% mutate(site = "Kalobeyei")
loturerei_daily_95th <- loturerei_daily_95th %>% mutate(site = "Loturerei")

all_sites_95head <- bind_rows(nadapal_daily_95th,kakemera_daily_95th,kalobeyei_daily_95th,loturerei_daily_95th)



##plot
all_sites_daily_95head_plot = ggplot(all_sites_95head, aes(x = date, y=p95_head/1000)) +
  geom_line()+
  labs(x = "",
       y = "95th percentile daily \n water column head (m)",
  ) +facet_wrap(~site, scales="free_y")+
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    axis.title.y = element_text(size = 12)
  )

all_sites_daily_95head_plot

ggsave("./plots/final-paper/all_sites_daily_95head_plot.png", plot = all_sites_daily_95head_plot, device = "png", dpi = 300,
       width = 6, height = 4, units = "in")

###TAMSAT monthly rainfall data - regionally applicable for all sites
tur_tamsat_monthly <- tur_tamsat_monthly %>% subset(date >= "2023-04-16" & date <= "2024-12-01")#"tur_tamsat_monthly" carried over from the "rainfall_seasonality_EAfrica.R" script
tur_tamsat_monthly <- tur_tamsat_monthly %>% mutate(year_mon = floor_date(date, "month"))

all_site_tamsat_monthly <- tur_tamsat_monthly %>% dplyr::select(year_mon,Kakemera,Loturerei,`Nadapal-IRC-1`,Kalobeyei)
colnames(all_site_tamsat_monthly)[4] <- "Nadapal-IRC"

all_site_tamsat_monthly <- all_site_tamsat_monthly %>%
  pivot_longer(cols = 2:5,
               names_to = "site", values_to = "rainfall")

#plot rainfall and corresponding groundwater heads
#all sites
all_sites_95head <- all_sites_95head %>% mutate(year_mon = floor_date(date, "month"))


rainfall_p95_head_plot <- ggplot() +
  # Daily groundwater head (line)
  geom_line(data = all_sites_95head, aes(x = date, y = p95_head/1000), color = "black",size=0.5) + #convert heads to m
  geom_point(data = all_sites_95head, aes(x = date, y = p95_head/1000), color = "black", size= 0.1)+
  
  # Monthly rainfall (bars)
  geom_col(data = all_site_tamsat_monthly, aes(x = year_mon, y = rainfall/10), 
           fill = "blue", alpha = 0.5) +
  
  facet_wrap(~ site, scales = "free_y", ncol = 2) +
  scale_y_continuous(
    name = "95th percentile daily \n water column head (m)", 
    sec.axis = sec_axis(~.*10, name = "Monthly rainfall (mm)"))+
  labs(x = "", y = "") +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("./plots/final-paper/f_v2/rainfall_p95_head_plot_v2.png", plot = rainfall_p95_head_plot, device = "png", dpi = 600,
       width = 6, height = 4, units = "in")

head(all_sites_95head)
