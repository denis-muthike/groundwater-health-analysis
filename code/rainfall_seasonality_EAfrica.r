##Rainfall anomalies analysis
##Prepared for the manuscript XXXX submitted to AGU Water Resources Research journal

#Load libraries, if not installed, use install.packages("") to install the following packages and load them into your R environment

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

setwd("/home/denis/projects/drip")

tur_gwl_sites <- st_read("./waterlevel/gwl-sites/tur_gwl_sites.shp")
tur_basins <- shapefile("./waterlevel/ancillary/turkana-riverbasins.shp")

##########################################################################################

###TAMSAT
# Load rainfall NetCDF (monthly rainfall in mm or similar units)
tamsat_daily_rain <- rast("./climate/tamsat/TAMSAT-v3p1-gapfilled-NKenya-2023-2024-daily-corrected.nc")
rain <- rast("./climate/tamsat/tamsat_v3p1_1983-2024-monthly.nc")
tam_dates <- seq(as.Date("1983-01-01"), as.Date("2024-12-01"), by="months")
names(rain)<-tam_dates
# Extract time dimension
# time_vals <- time(rain)
# dates <- as.Date(time_vals)
# names(rain) <- format(dates, "%Y-%m")

# Seasons definition
season_def <- list(
  MAM  = c(3,4,5),       # long rains
  JJAS = c(6,7,8,9),     # june-sept
  OND  = c(10,11,12)     # short rains
)

# Subset climatology period (1983–2012)
clim_period <- tam_dates >= as.Date("1983-01-01") & tam_dates <= as.Date("2024-12-31")
rain_clim <- rain[[clim_period]]
dates_clim <- tam_dates[tam_dates >= "1983-01-01" & tam_dates <= "2024-12-31"]

# Function to compute seasonal means for a given period
compute_seasonal <- function(r, dates, season_months) {
  df <- tibble(date = dates,
               year = year(dates),
               month = month(dates))
  
  # Identify seasonal groups
  df <- df %>%
    dplyr::filter(month %in% season_months) %>%
    group_by(year) %>%
    summarise(idx = list(which(month(dates) %in% season_months & year(dates)==year)))
  
  # Stack seasonal sums
  out_list <- lapply(seq_len(nrow(df)), function(i) {
    sum(r[[df$idx[[i]]]]) # seasonal total/mean (here SUM monthly within season, can use MEAN if needed)
  })
  
  s <- rast(out_list)
  names(s) <- df$year
  return(s)
}

# --- Compute climatology for each season (1983–2024) ---
climatologies <- lapply(season_def, function(m) compute_seasonal(rain_clim, dates_clim, m))

# Average climatology (long-term mean per pixel)
clim_mean <- lapply(climatologies, function(r) mean(r, na.rm=TRUE))

# --- Compute anomalies for 2000–2024 ---
analysis_period <- tam_dates >= as.Date("2000-01-01")
rain_analysis <- rain[[analysis_period]]
dates_analysis <- tam_dates[tam_dates >= "2000-01-01"]

seasonal_analysis <- lapply(names(season_def), function(season) {
  # seasonal totals for each year
  seas_stack <- compute_seasonal(rain_analysis, dates_analysis, season_def[[season]])
  
  # match climatology mean for this season
  clim <- clim_mean[[season]]
  
  # anomalies as % of long-term average
  anomaly <- (seas_stack / clim) * 100
  names(anomaly) <- names(seas_stack)
  
  return(anomaly)
})
names(seasonal_analysis) <- names(season_def)

# Example: plot anomaly for OND in 2019
plot(seasonal_analysis$OND[["2019"]], main="OND 2019 Anomaly (% of 1983–2024 mean)")


#############################################################
####Regional rainfall anomalies - using basins to extract
###############################################################

subset_basins <- st_as_sf(tur_basins)
subset_basins <- subset_basins %>% dplyr::filter(LBASIN_ID %in% c(1,3,35,47))#retain only desired features
subset_basins <- subset_basins %>% summarise() #dissolve/union into one feature

Seasonal_anomalies_OND <- seasonal_analysis[["OND"]]
Seasonal_anomalies_MAM <- seasonal_analysis[["MAM"]]
Seasonal_anomalies_JJAS <- seasonal_analysis[["JJAS"]]

tur_basins_OND_anom <- raster::extract(Seasonal_anomalies_OND, subset_basins, mean, na.rm=T)

tur_basins_OND_anom_long <- tur_basins_OND_anom %>%
  pivot_longer(
    cols = -ID,              # keep Site as is
    names_to = "Year",   # column names go here
    values_to = "ppt_anomaly"
  )

head(tur_basins_OND_anom_long)

tur_basins_OND_anom_long <- tur_basins_OND_anom_long %>%
  mutate(perc_anom_dev = ppt_anomaly - 100)  # deviation from 100%

ggplot(tur_basins_OND_anom_long, aes(x = Year, y = perc_anom_dev)) +
  geom_col(position = "dodge") +               # side-by-side bars
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # baseline = 100%
  labs(title = "OND Rainfall Anomalies (relative to 1983–2024 mean)",
       y = "Deviation from long-term mean (%)",
       x = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))

###MAM
tur_basins_MAM_anom <- raster::extract(Seasonal_anomalies_MAM, subset_basins, mean, na.rm=T)

tur_basins_MAM_anom_long <- tur_basins_MAM_anom %>%
  pivot_longer(
    cols = -ID,              # keep Site as is
    names_to = "Year",   # column names (MAM2000 etc.) go here
    values_to = "ppt_anomaly"
  )

head(tur_basins_MAM_anom_long)

tur_basins_MAM_anom_long <- tur_basins_MAM_anom_long %>%
  mutate(perc_anom_dev = ppt_anomaly - 100)  # deviation from 100%

ggplot(tur_basins_MAM_anom_long, aes(x = Year, y = perc_anom_dev)) +
  geom_col(position = "dodge") +               # side-by-side bars
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # baseline = 100%
  labs(title = "MAM Rainfall Anomalies (relative to 1983–2024 mean)",
       y = "Percent anomaly",
       x = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5))

###JJAS
tur_basins_JJAS_anom <- raster::extract(Seasonal_anomalies_JJAS, subset_basins, mean, na.rm=T)

tur_basins_JJAS_anom_long <- tur_basins_JJAS_anom %>%
  pivot_longer(
    cols = -ID,              # keep Site as is
    names_to = "Year",   # column names (JJAS2000 etc.) go here
    values_to = "ppt_anomaly"
  )

head(tur_basins_JJAS_anom_long)

tur_basins_JJAS_anom_long <- tur_basins_JJAS_anom_long %>%
  mutate(perc_anom_dev = ppt_anomaly - 100)  # deviation from 100%

ggplot(tur_basins_JJAS_anom_long, aes(x = Year, y = perc_anom_dev)) +
  geom_col(position = "dodge") +               # side-by-side bars
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # baseline = 100%
  labs(title = "JJAS Rainfall Anomalies (relative to 1983–2024 mean)",
       y = "Anomaly (%)",
       x = "Year") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

###Combine
tur_basins_MAM_anom_long <- tur_basins_MAM_anom_long %>% mutate(Season = "MAM")
tur_basins_JJAS_anom_long <- tur_basins_JJAS_anom_long %>% mutate(Season = "JJAS")
tur_basins_OND_anom_long <- tur_basins_OND_anom_long %>% mutate(Season = "OND")

all_seasonal_rainfall_anomalies <- rbind(tur_basins_MAM_anom_long,tur_basins_OND_anom_long)

#plot
all_seasonal_rainfall_anom_plot = ggplot(all_seasonal_rainfall_anomalies %>% dplyr::filter(Year >= 2002 & Year <= 2024), aes(x = Year, y = perc_anom_dev, fill = Season, alpha = Season)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("MAM" = "darkgreen", 
                               "OND" = "blue")) +
  scale_alpha_manual(values = c("MAM" = 1, 
                                "OND" = 1)) +
  labs(x = "",
       y = "Seasonal rainfall anomaly (%) \n based on the 1983-2024 mean"
  ) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.background = element_rect(
          fill = "white",        # legend background color
          color = "black",       # outline color
          linewidth = 0.1),       # thickness of border
        legend.position = c(0.9, 0.85),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12))

ggsave("./plots/final-paper/f_v2/tamsat-all-sites-rainfall-anomalies_v3.png", plot = all_sites_rainfall_anomalies,
       device = png,
       dpi = 500,
       width = 6, height = 4, units = "in")           # No scaling



####Gridded maps for paper
########################################################################################3
###many pluvial and drought seasons

########################################################################################3



mult_OND <- Seasonal_anomalies_OND[[c("2005","2010","2011","2015","2019")]]
mult_MAM <- Seasonal_anomalies_MAM[[c("2010","2011","2018","2020","2023")]]

mult_OND_anom <- terra::crop(mult_OND,tur_basins)
mult_MAM_anom <- terra::crop(mult_MAM,tur_basins)

names(mult_OND_anom) <- c("OND_2005","OND_2010", "OND_2011", "OND_2015","OND_2019")
names(mult_MAM_anom) <- c("MAM_2010","MAM_2011", "MAM_2018", "MAM_2020","MAM_2023")

stack_anom <- c(mult_MAM_anom,mult_OND_anom)

myTheme_1 <- rasterTheme(region = colorRampPalette(c("red","yellow","green", "darkgreen","lightblue","blue"))(200))
myTheme_2 <- rasterTheme(region = colorRampPalette(c("brown", "yellow","blue"))(100))

myTheme_1_with_border <- modifyList(myTheme_1, list(
  axis.line = list(col = "black", lwd = 1.5)
))

breaks <- c(0, 50, 100, 150, max(values(stack_anom), na.rm = TRUE))
colors <- c("tan","yellow", "green","lightblue", "darkblue") 

###resample for smoothed viz
# Build a template grid at the resolution you want
tgt <- rast(stack_anom)
res(tgt) <- c(0.001, 0.001)   # lon, lat degrees
# Optional: align origin
origins <- origin(stack_anom); origin(tgt) <- origins

# Resample with bilinear (continuous) or near (categorical)
x_target <- resample(stack_anom, tgt, method = "bilinear")

# x_target <- raster::stack(x_target)

rain_plot = levelplot(x_target,main = "Seasonal Rainfall Anomalies",margin=FALSE,par.settings=myTheme_1_with_border,colorkey=list(space="bottom",title=list(label = "% of avg.",side = "right",rot = 0),
                                                                                   height = 0.8,par.strip.text = list(cex = 0.5)),xlab = NULL, ylab=NULL)

rain_plot_final = rain_plot+latticeExtra::layer(sp.lines(tur_basins,lwd=1,col="black"))

png("./plots/final-paper/f_v2/tamsat-seasonal-anomaly-maps-v4.png", width = 3000, height = 2000, res = 350)
rain_plot_final
dev.off()

pdf("./plots/final-paper/tamsat-seasonal-anomaly-maps.pdf", width = 10, height = 7)
rain_plot_final
dev.off()


#####Site-level rainfall-waterlevel correlations
tamsat_mon <- rast("./climate/tamsat/tamsat_v3p1_1983-2024-monthly.nc")

tur_tamsat_daily <- raster::extract(tamsat_daily_rain, tur_gwl_sites, mean, na.rm=T)

tamsat_daily_dates <- seq(as.Date("2023-01-01"), as.Date("2024-12-31"), by="days")
tamsat_monthly_dates<- seq(as.Date("1983-01-01"), as.Date("2024-12-31"), by="months")

tur_tamsat_daily <- as.data.frame(t(tur_tamsat_daily))
tur_tamsat_daily <- tur_tamsat_daily[-1,]
tur_tamsat_daily <- rownames_to_column(tur_tamsat_daily, var="date")
tur_tamsat_daily[1]<-tamsat_daily_dates
colnames(tur_tamsat_daily)[2:7] <-tur_gwl_sites$Name

#monthly
tur_tamsat_monthly <- raster::extract(tamsat_mon, tur_gwl_sites, mean, na.rm=T)
tur_tamsat_monthly <- as.data.frame(t(tur_tamsat_monthly))
tur_tamsat_monthly <- tur_tamsat_monthly[-1,]
tur_tamsat_monthly <- rownames_to_column(tur_tamsat_monthly, var="date")
tur_tamsat_monthly[1]<-tamsat_monthly_dates
colnames(tur_tamsat_monthly)[2:7] <-tur_gwl_sites$Name

mean_mon_rainfall <- tur_tamsat_monthly %>% transmute(date, mean_rainfall = rowMeans(dplyr::select(., -date), na.rm = TRUE))

#Loturerei
loturerei_tamsat_rf <- tur_tamsat_daily %>% dplyr::select(date,Loturerei)
colnames(loturerei_tamsat_rf)[2]<-"rf_1d"
loturerei_tamsat_rf$rf_1d <- as.numeric(loturerei_tamsat_rf$rf_1d)

loturerei_tamsat_wk_rf <- loturerei_tamsat_rf %>% dplyr::select(date,rf_1d)
loturerei_week_rf <- loturerei_tamsat_wk_rf %>% 
  mutate(week = isoweek(date),
         year = year(date)) %>% group_by(year,week) %>% 
  summarise(total_rf = sum(rf_1d, na.rm = T),
            mean_rf = mean(rf_1d, na.rm = T),
            .groups = "drop")

#cumulative rainfall
loturerei_tamsat_rf <- setDF(loturerei_tamsat_rf) #convert to data.table format
cumulative_windows <- c(10,15,30,60,90,120)

#rolling sums
for (k in cumulative_windows){
  loturerei_tamsat_rf[[paste0("rf_",k)]] <- frollsum(loturerei_tamsat_rf$rf_1d, k, fill = NA,align = "right")
}

loturerei_95head <- loturerei_daily_95th #pulled from "groundwater-head-trend-analysis.R"
loturerei_rainfall_gwlhead <- merge(loturerei_tamsat_rf,loturerei_95head, by="date")

#correlations
loturerei_correlations <- loturerei_rainfall_gwlhead %>% dplyr::select(p95_head,starts_with("rf"))

loturerei_corr_mat <- cor(loturerei_correlations, use = "pairwise.complete.obs")
print(loturerei_corr_mat)

vars <- names(loturerei_correlations)
n <- length(vars)
corr_mat <- matrix(NA_real_, n, n, dimnames = list(vars, vars))
p_mat    <- matrix(NA_real_, n, n, dimnames = list(vars, vars))

for (i in seq_len(n)) {
  for (j in i:n) {
    x <- loturerei_correlations[[i]]
    y <- loturerei_correlations[[j]]
    ok <- complete.cases(x, y)
    ct <- cor.test(x[ok], y[ok], method = "pearson")  # or "spearman"/"kendall"
    corr_mat[i, j] <- corr_mat[j, i] <- ct$estimate
    p_mat[i, j]    <- p_mat[j, i]    <- ct$p.value
  }
}

corr_mat
p_mat

#Kalobeyei
kalobeyei_tamsat_rf <- tur_tamsat_daily %>% dplyr::select(date,Kalobeyei)
colnames(kalobeyei_tamsat_rf)[2]<-"rf_1d"
kalobeyei_tamsat_rf$rf_1d <- as.numeric(kalobeyei_tamsat_rf$rf_1d)

#cumulative rainfall
kalobeyei_tamsat_rf <- setDF(kalobeyei_tamsat_rf) #convert to data.table format
cumulative_windows <- c(10,15,30,60,90,120)

#rolling sums
for (k in cumulative_windows){
  kalobeyei_tamsat_rf[[paste0("rf_",k)]] <- frollsum(kalobeyei_tamsat_rf$rf_1d, k, fill = NA,align = "right")
}

kalobeyei_95head <- kalobeyei_daily_95th #pulled from "groundwater-head-trend-analysis.R"
kalobeyei_rainfall_gwlhead <- merge(kalobeyei_tamsat_rf,kalobeyei_95head, by="date")

#correlations
kalobeyei_correlations <- kalobeyei_rainfall_gwlhead %>% dplyr::select(p95_head,starts_with("rf"))

kalobeyei_corr_mat <- cor(kalobeyei_correlations, use = "pairwise.complete.obs")
print(kalobeyei_corr_mat)

vars <- names(kalobeyei_correlations)
n <- length(vars)
corr_mat <- matrix(NA_real_, n, n, dimnames = list(vars, vars))
p_mat    <- matrix(NA_real_, n, n, dimnames = list(vars, vars))

for (i in seq_len(n)) {
  for (j in i:n) {
    x <- kalobeyei_correlations[[i]]
    y <- kalobeyei_correlations[[j]]
    ok <- complete.cases(x, y)
    ct <- cor.test(x[ok], y[ok], method = "pearson")  # or "spearman"/"kendall"
    corr_mat[i, j] <- corr_mat[j, i] <- ct$estimate
    p_mat[i, j]    <- p_mat[j, i]    <- ct$p.value
  }
}

corr_mat
p_mat

#Nadapal-IRC
nadapal_tamsat_rf <- tur_tamsat_daily %>% dplyr::select(date,5) #index "5" is Nadapal because of the -1 on the name label
colnames(nadapal_tamsat_rf)[2]<-"rf_1d"
nadapal_tamsat_rf$rf_1d <- as.numeric(nadapal_tamsat_rf$rf_1d)

#cumulative rainfall
nadapal_tamsat_rf <- setDF(nadapal_tamsat_rf) #convert to data.table format
cumulative_windows <- c(10,15,30,60,90,120)

#rolling sums
for (k in cumulative_windows){
  nadapal_tamsat_rf[[paste0("rf_",k)]] <- frollsum(nadapal_tamsat_rf$rf_1d, k, fill = NA,align = "right")
}

nadapal_95head <- nadapal_daily_95th #pulled from "groundwater-head-trend-analysis.R"
nadapal_rainfall_gwlhead <- merge(nadapal_tamsat_rf,nadapal_95head, by="date")

#correlations
nadapal_correlations <- nadapal_rainfall_gwlhead %>% dplyr::select(p95_head,starts_with("rf"))

nadapal_corr_mat <- cor(nadapal_correlations, use = "pairwise.complete.obs")
print(nadapal_corr_mat)

vars <- names(nadapal_correlations)
n <- length(vars)
corr_mat <- matrix(NA_real_, n, n, dimnames = list(vars, vars))
p_mat    <- matrix(NA_real_, n, n, dimnames = list(vars, vars))

for (i in seq_len(n)) {
  for (j in i:n) {
    x <- nadapal_correlations[[i]]
    y <- nadapal_correlations[[j]]
    ok <- complete.cases(x, y)
    ct <- cor.test(x[ok], y[ok], method = "pearson")  # or "spearman"/"kendall"
    corr_mat[i, j] <- corr_mat[j, i] <- ct$estimate
    p_mat[i, j]    <- p_mat[j, i]    <- ct$p.value
  }
}

corr_mat
p_mat



#Kakemera
kakemera_tamsat_rf <- tur_tamsat_daily %>% dplyr::select(date,Kakemera) 
colnames(kakemera_tamsat_rf)[2]<-"rf_1d"
kakemera_tamsat_rf$rf_1d <- as.numeric(kakemera_tamsat_rf$rf_1d)

#cumulative rainfall
kakemera_tamsat_rf <- setDF(kakemera_tamsat_rf) #convert to data.table format
cumulative_windows <- c(10,15,30,60,90,120)

#rolling sums
for (k in cumulative_windows){
  kakemera_tamsat_rf[[paste0("rf_",k)]] <- frollsum(kakemera_tamsat_rf$rf_1d, k, fill = NA,align = "right")
}

kakemera_95head <- kakemera_daily_95th #pulled from "groundwater-head-trend-analysis.R"
kakemera_rainfall_gwlhead <- merge(kakemera_tamsat_rf,kakemera_95head, by="date")

#correlations
kakemera_correlations <- kakemera_rainfall_gwlhead %>% dplyr::select(p95_head,starts_with("rf"))

kakemera_corr_mat <- cor(kakemera_correlations, use = "pairwise.complete.obs")
print(kakemera_corr_mat)

vars <- names(kakemera_correlations)
n <- length(vars)
corr_mat <- matrix(NA_real_, n, n, dimnames = list(vars, vars))
p_mat    <- matrix(NA_real_, n, n, dimnames = list(vars, vars))

for (i in seq_len(n)) {
  for (j in i:n) {
    x <- kakemera_correlations[[i]]
    y <- kakemera_correlations[[j]]
    ok <- complete.cases(x, y)
    ct <- cor.test(x[ok], y[ok], method = "pearson")  # or "spearman"/"kendall"
    corr_mat[i, j] <- corr_mat[j, i] <- ct$estimate
    p_mat[i, j]    <- p_mat[j, i]    <- ct$p.value
  }
}

corr_mat
p_mat



###combine daily head and rainfall - only the 4 selected sites

loturerei_rainfall_gwlhead <- loturerei_rainfall_gwlhead %>% mutate(site = "Loturerei")
kakemera_rainfall_gwlhead <- kakemera_rainfall_gwlhead %>% mutate(site = "Kakemera")
kalobeyei_rainfall_gwlhead <- kalobeyei_rainfall_gwlhead %>% mutate(site = "Kalobeyei")
nadapal_rainfall_gwlhead <- nadapal_rainfall_gwlhead %>% mutate(site = "Nadapal-IRC")

all_sites_rainfall_gwlhead <- bind_rows(nadapal_rainfall_gwlhead,kalobeyei_rainfall_gwlhead,kakemera_rainfall_gwlhead,loturerei_rainfall_gwlhead)

all_sites_rainfall_gwlhead <- all_sites_rainfall_gwlhead %>% dplyr::select(date,rf_1d,p95_head,site)



#####Terrestrial Water Availability
##Source: GRACE/GRACE-FO: JPL & CSR; https://grace.jpl.nasa.gov/data-analysis-tool

tws_anomaly <- turkana_grace_tws_df %>% dplyr::select(year_mon,mean_tws)#carried over from the "groundwater-analysis.R" script
colnames(tws_anomaly)[1]<-"date"

#merge
rain_tws_anom <- merge(mon_rain_anom,tws_anomaly, by="date")

rain_tws_anom <- rain_tws_anom %>% mutate(year = year(date)) %>% group_by(year) %>% 
                                        summarise(mean_tws = mean(jpl_lwe, na.rm=T),
                                                  mean_rf_pct = mean(anom_pct, na.rm=T))

tws_anomaly_plot = ggplot(tws_anomaly, aes(x=datetime, y=jpl_lwe))+geom_area(color="blue",fill="blue")+
  geom_hline(yintercept = 0, linetype = "dashed")+labs(x="",y="Water Equivalent \n Thickness (cm)")+theme_minimal(base_size = 12)


ggplot(rain_tws_anom, aes(x = date)) +
  # TWS anomaly (left axis)
  geom_line(aes(y = mean_tws, color = "TWS anomaly"), size = 1) +
  # Rain anomaly scaled to match the second axis
  geom_line(aes(y = anom_pct, color = "Rainfall anomaly"), size = 1) +

  scale_y_continuous(
    name = "TWS anomaly (cm)",
    sec.axis = sec_axis(~ ., name = "Rainfall anomaly (mm)") # adjust if scaling applied
  ) +
  scale_color_manual(
    name = "",
    values = c("TWS anomaly" = "blue", "Rainfall anomaly" = "brown")
  ) +
  labs(
    x = "Date",
    title = "TWS and Rainfall Anomalies"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    axis.title.y.left = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "brown")
  )

tws_anomaly_plot = ggplot(tws_anomaly, aes(x = datetime)) +
                  geom_ribbon(aes(ymin = pmin(jpl_lwe, 0), ymax = 0,
                  fill = "Below Normal"), alpha = 1,show.legend = FALSE) +
                  geom_ribbon(aes(ymin = 0, ymax = pmax(jpl_lwe, 0),
                  fill = "Above Normal"), alpha = 1,show.legend = FALSE) +  
                  geom_smooth(aes(y=jpl_lwe),method = "loess", se = FALSE,
                  color = "black", linewidth = 1.3, span = 0.3) +
                  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
                  scale_fill_manual(values = c("Above Normal" = "blue", "Below Normal" = "brown"),
                  name = "Anomaly") +
                  labs(x="",y="GRACE water equivalent \n thickness anomaly (cm)") +
                  theme_bw()+theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 0.5),
                                   axis.text.y = element_text(size = 10),
                                   axis.title.y = element_text(size = 12))


ggsave("D:/drip/plots/final-paper/f_v2/grace-tws_anomaly_plot_3.png", plot = tws_anomaly_plot, device = "png", dpi = 350,
       width = 6, height = 4, units = "in")

###combine plots for rainfall and tws anomalies
combined_plot = grid.arrange(all_seasonal_rainfall_anom_plot,tws_anomaly_plot, nrow=2)

ggsave("D:/drip/plots/final-paper/f_v2/rainfall-and-tws-anomalies-v1.png", plot = combined_plot, device = "png", dpi = 300,
       width = 6, height = 6, units = "in")


##################################
###Frequency of wet/dry seasons
####SPI-based analysis
####################################################
# Load monthly TAMSAT data
rain <- rast("./climate/tamsat/tamsat_v3p1_1983-2024-monthly.nc")
rain_clip <- terra::crop(rain,tur_basins)

# Get time information and create year-month dataframe
dates <- tam_dates
years <- year(dates)
months <- month(dates)

cat("Data range:", as.character(min(dates)), "to", as.character(max(dates)), "\n")
cat("Total months:", nlyr(rain_clip), "\n")

# Calculate SPI for main rainy seasons
calculate_key_seasons_spi <- function(monthly_rast, dates) {
  years <- year(dates)
  months <- month(dates)
  
  key_seasons <- list(
    MAM = 3:5,   # Long Rains - proper 3-month accumulation
    JJAS = 6:9,  # Long Dry - 4-month accumulation (needs special handling)
    OND = 10:12  # Short Rains - proper 3-month accumulation
  )
  
  seasonal_data <- list()
  season_info <- data.frame()
  
  for(season_name in names(key_seasons)) {
    season_months <- key_seasons[[season_name]]
    cat("Processing", season_name, "with", length(season_months), "month accumulation...\n")
    
    seasonal_layers <- list()
    
    for(yr in unique(years)) {
      month_mask <- years == yr & months %in% season_months
      
      if(sum(month_mask) == length(season_months)) {
        # Calculate seasonal total
        season_total <- sum(monthly_rast[[which(month_mask)]])
        seasonal_layers[[paste0(yr, "_", season_name)]] <- season_total
        
        season_info <- rbind(season_info,
                             data.frame(Year = yr,
                                        Season = season_name,
                                        LayerName = paste0(yr, "_", season_name)))
      }
    }
    
    if(length(seasonal_layers) > 10) {
      seasonal_rast <- rast(seasonal_layers)
      seasonal_data[[season_name]] <- seasonal_rast
    }
  }
  
  return(list(seasonal_data = seasonal_data, season_info = season_info))
}

# Calculate key seasons with proper accumulation

seasonal_result <- calculate_key_seasons_spi(rain_clip, dates)
seasonal_data <- seasonal_result$seasonal_data
season_info <- seasonal_result$season_info

# CALCULATE SPI FOR EACH SEASON

calculate_seasonal_zscore_from_totals <- function(seasonal_data, season_info) {
  zscore_results <- list()
  
  for (season_name in names(seasonal_data)) {
    cat("Calculating zscore for", season_name, "...\n")
    season_rast <- seasonal_data[[season_name]]
    
    # Calculate zscore
    zscore_rast <- app(season_rast, function(x) {
      if (sum(!is.na(x)) >= 10) {
        (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
      } else {
        NA
      }
    })
    
    names(zscore_rast) <- names(season_rast)
    zscore_results[[season_name]] <- zscore_rast
    
    # Print zscore statistics for quality control
    zscore_values <- values(zscore_rast)
    cat("  zscore range:", round(min(zscore_values, na.rm = TRUE), 2), "to",
        round(max(zscore_values, na.rm = TRUE), 2), "\n")
    cat("  Mean zscore:", round(mean(zscore_values, na.rm = TRUE), 2), "\n")
    cat("  SD of zscore:", round(sd(zscore_values, na.rm = TRUE), 2), "\n")
  }
  
  return(zscore_results)
}
# Calculate zscore for each season
zscore_results <- calculate_seasonal_zscore_from_totals(seasonal_data, season_info)

# GROUPED CATEGORIES ANALYSIS

analyze_grouped_zscore_categories <- function(zscore_results, season_info) {
  extreme_counts <- list()
  temporal_data <- data.frame()
  
  for(season_name in names(zscore_results)) {
    zscore_rast <- zscore_results[[season_name]]
    
    # Get corresponding years
    season_years <- season_info$Year[season_info$Season == season_name]
    
    # Calculate extreme counts for each category
    # DROUGHTS: zscore ≤ -0.8 (D1-D4)
    drought_count <- app(zscore_rast, function(x) sum(x <= -0.8, na.rm = TRUE))
    
    # FLOODS: zscore ≥ 0.8 (W2-W5)
    flood_count <- app(zscore_rast, function(x) sum(x >= 0.8, na.rm = TRUE))
    
    # SEVERE DROUGHTS: zscore ≤ -1.3 (D2-D4)
    severe_drought_count <- app(zscore_rast, function(x) sum(x <= -1.3, na.rm = TRUE))
    
    # SEVERE FLOODS: zscore ≥ 1.3 (W3-W5)
    severe_flood_count <- app(zscore_rast, function(x) sum(x >= 1.3, na.rm = TRUE))
    
    extreme_counts[[paste0(season_name, "_drought")]] <- drought_count
    extreme_counts[[paste0(season_name, "_flood")]] <- flood_count
    extreme_counts[[paste0(season_name, "_severe_drought")]] <- severe_drought_count
    extreme_counts[[paste0(season_name, "_severe_flood")]] <- severe_flood_count
    
    # Calculate temporal evolution
    for(i in 1:nlyr(zscore_rast)) {
      year <- season_years[i]
      current_zscore <- zscore_rast[[i]]
      
      # Calculate area affected
      drought_area <- global(current_zscore <= -0.8, mean, na.rm = TRUE)[1,1] * 100
      flood_area <- global(current_zscore >= 0.8, mean, na.rm = TRUE)[1,1] * 100
      severe_drought_area <- global(current_zscore <= -1.3, mean, na.rm = TRUE)[1,1] * 100
      severe_flood_area <- global(current_zscore >= 1.3, mean, na.rm = TRUE)[1,1] * 100
      
      temporal_data <- rbind(temporal_data,
                             data.frame(Year = year,
                                        Season = season_name,
                                        Drought_Area = drought_area,
                                        Flood_Area = flood_area,
                                        Severe_Drought_Area = severe_drought_area,
                                        Severe_Flood_Area = severe_flood_area))
    }
    
    cat("Processed", season_name, "-", nlyr(zscore_rast), "years\n")
  }
  
  return(list(extreme_counts = rast(extreme_counts),
              temporal_data = temporal_data))
}

# Perform grouped analysis
grouped_analysis <- analyze_grouped_zscore_categories(zscore_results, season_info)
extreme_counts <- grouped_analysis$extreme_counts
temporal_data <- grouped_analysis$temporal_data

# VISUALIZE RESULTS

##subset to MAM and OND only
extreme_counts_subset <- extreme_counts[[c("MAM_drought","MAM_flood","OND_drought","OND_flood")]]

##drought
drought_subset <- extreme_counts_subset[[c("MAM_drought","OND_drought")]]
names(drought_subset)<-c("MAM","OND")
drought_subset <- 100*(drought_subset/42)
flood_subset <- extreme_counts_subset[[c("MAM_flood","OND_flood")]]
names(flood_subset)<-c("MAM","OND")
flood_subset <- 100*(flood_subset/42)

# Build a template grid at a higher spatial resolution
tgt_drg <- rast(drought_subset)
res(tgt_drg) <- c(0.001, 0.001)   # lon, lat degrees
# Optional: align origin
origins <- origin(drought_subset); origin(tgt_drg) <- origins

# Resample with bilinear (continuous) or near (categorical)
x_target_drgt <- resample(drought_subset, tgt_drg, method = "bilinear")

####flood
# Build a template grid at the resolution you want
tgt_fld <- rast(flood_subset)
res(tgt_fld) <- c(0.001, 0.001)   # lon, lat degrees
# Optional: align origin
origins <- origin(flood_subset); origin(tgt_fld) <- origins

# Resample with bilinear (continuous) or near (categorical)
x_target_fld <- resample(flood_subset, tgt_fld, method = "bilinear")

cols_flood <- hcl.colors(50, "Blues", rev = TRUE)
cols_drought <- hcl.colors(20, "Reds", rev = TRUE)

##alternative colors
r_pos <- terra::app(x_target_drgt, function(x) pmax(x, 0))
rmax   <- max(values(r_pos), na.rm = TRUE)
at     <- c(5, 10, 15, 20, 25, 30, rmax)   # color intervals

drought_freq_plot = levelplot(x_target_drgt,main = "Frequency of dry events (% of climatology)",margin=FALSE,at = at, col.regions=cols_drought,colorkey=list(space="bottom",at = at, title=list(label = "Frequency (%)", side = "bottom", rot = 0),
                                                                                                                                    height = 0.8,par.strip.text = list(cex = 0.5)),xlab = NULL, ylab=NULL)

drought_freq_plot_final = drought_freq_plot+latticeExtra::layer(sp.lines(tur_basins,lwd=1,col="black"))

png("./plots/final-paper/f_v2/drought_freq_plot_final.png", width = 3000, height = 2000, res = 350)
drought_freq_plot_final
dev.off()


##flood
r_pos <- terra::app(x_target_fld, function(x) pmax(x, 0))

# For a proper legend tick at the last bin, cap Inf at the raster max:
rmax_fld   <- max(values(r_pos), na.rm = TRUE)
at_fld     <- c(0, 5, 10, 15, 20, 25, 30, rmax_fld)   # color intervals


flood_freq_plot = levelplot(x_target_fld,main = "Frequency of wet events (% of climatology)",margin=FALSE,col.regions=cols_flood, at = at_fld, colorkey=list(space="bottom",at = at_fld,title=list(label = "Frequency (%)", side = "bottom", rot = 0),
                                                                                                                            height = 0.8,par.strip.text = list(cex = 0.5)),xlab = NULL, ylab=NULL)

flood_freq_plot_final = flood_freq_plot+latticeExtra::layer(sp.lines(tur_basins,lwd=1,col="black"))

png("./plots/final-paper/f_v2/flood_freq_plot_final.png", width = 3000, height = 2000, res = 350)
flood_freq_plot_final
dev.off()


############################
#CDO based-seasonal-rainfall-anomaly-analysis

#####
MAM_seas_perc_anom <- rast("./climate/MAM-percent-of-climatology-1983-2024.nc")
OND_seas_perc_anom <- rast("./climate/OND-percent-of-climatology-1983-2024.nc")

names(MAM_seas_perc_anom)<-1983:2024
names(OND_seas_perc_anom)<-1983:2024

yrs <- as.integer(trimws(names(MAM_seas_perc_anom)))  # same for OND

MAM_select <- 2000:2024#c(2010, 2011, 2018, 2020, 2023)
OND_select <- c(2005, 2010, 2011, 2015, 2019)

MAM_select_anom <- MAM_seas_perc_anom[[ yrs %in% MAM_select ]]

OND_select_anom <- OND_seas_perc_anom[[ yrs %in% OND_select ]]

tgt <- rast(MAM_select_anom)
res(tgt) <- c(0.001, 0.001)   # lon, lat degrees
# Optional: align origin
origins <- origin(MAM_select_anom); origin(tgt) <- origins

ond_tgt <- rast(OND_select_anom)
res(ond_tgt) <- c(0.001, 0.001)   # lon, lat degrees
# Optional: align origin
ond_origins <- origin(OND_select_anom); origin(ond_tgt) <- origins

# Resample with bilinear (continuous) or near (categorical)
mam_target <- resample(MAM_select_anom, tgt, method = "bilinear")

names(mam_target) <- c("MAM_2010","MAM_2011", "MAM_2018","MAM_2020","MAM_2023")

mam_target <- terra::crop(mam_target, tur_basins)

mam_target_stk <- stack(mam_target)

##OND
ond_target <- resample(OND_select_anom, ond_tgt, method = "bilinear")

names(ond_target) <- c("OND_2005","OND_2010", "OND_2011","OND_2015","OND_2019")

ond_target <- terra::crop(ond_target, tur_basins)

ond_target_stk <- stack(ond_target)

##combined stack

multi_seas_stk <- stack(mam_target_stk,ond_target_stk)


########################################################################################3
###many pluvial and drought seasons

########################################################################################3

myTheme_1 <- rasterTheme(region = colorRampPalette(c("#8C510A","#D8B365","lightgreen","green", "darkgreen","lightblue","blue"))(200))

myTheme_1_with_border <- modifyList(myTheme_1, list(
  axis.line = list(col = "black", lwd = 1.5)
))

##
rain_pal <- colorRampPalette(
  c("#8C510A", "#D8B365", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#5AB4AC", "#01665E")
)(50)

myTheme_rain_anom <- rasterTheme(region = rain_pal)

myTheme_2_with_border <- modifyList( myTheme_rain_anom, list(
  axis.line = list(col = "black", lwd = 1.5)
))

new_rain_plot = levelplot(multi_seas_stk,main = "Seasonal Rainfall Anomalies",margin=FALSE,par.settings=myTheme_1_with_border,colorkey=list(space="bottom",title=list(label = "% of avg.",side = "right",rot = 0),
                                                                                                                                  height = 0.8,par.strip.text = list(cex = 0.5)),xlab = NULL, ylab=NULL)

rain_plot_final = new_rain_plot+latticeExtra::layer(sp.lines(tur_basins,lwd=1,col="black"))

png("./plots/final-paper/f_v2/tamsat-seasonal-anomaly-maps-final.png", width = 3000, height = 2000, res = 350)
rain_plot_final
dev.off()

pdf("./plots/final-paper/tamsat-seasonal-anomaly-maps.pdf", width = 10, height = 7)
rain_plot_final
dev.off()
