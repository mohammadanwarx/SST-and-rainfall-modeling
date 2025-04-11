

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
install.packages(c("tidyverse", "gifski", "terra", "gganimate"))
library(terra)
library(ggplot2)

# Load datasets
sst_file <- "D:/Geo x/SST and rainfall/sst.mon.anom.nc"
sst_raster <- rast(sst_file)

rain_file <- "D:/Geo x/SST and rainfall/sudan_rain_1981_2021.nc"
rain_raster <- rast(rain_file)

# Process Niño 3.4 Region
unzip("D:/Geo x/SST and rainfall/Nini3.ZIP", exdir = "D:/Geo x/SST and rainfall/Nini3")
shp_file <- list.files("D:/Geo x/SST and rainfall/Nini3", pattern = "\\.shp$", full.names = TRUE)
roi <- vect(shp_file)

# Adjust SST longitude from 0-360 to -180-180 and align CRS
sst_rotated <- rotate(sst_raster)
if (crs(roi) != crs(sst_rotated)) {
  roi <- project(roi, crs(sst_rotated))
}
sst_clipped <- crop(sst_rotated, roi)
sst_clipped <- mask(sst_clipped, roi)

# Process Khartoum Rainfall
unzip("D:/Geo x/SST and rainfall/Khartoum.ZIP", exdir = "D:/Geo x/SST and rainfall/Khartoum_shp")
khartoum_shp <- vect("D:/Geo x/SST and rainfall/Khartoum_shp/Khartoum.shp")

# Align CRS for rainfall data
if (crs(khartoum_shp) != crs(rain_raster)) {
  khartoum_shp <- project(khartoum_shp, crs(rain_raster))
}
rain_clipped <- crop(rain_raster, khartoum_shp)
rain_clipped <- mask(rain_clipped, khartoum_shp)

# Align dates between SST and rainfall
rain_dates <- time(rain_clipped)
sst_clipped <- sst_clipped[[which(time(sst_clipped) %in% rain_dates)]]

# Verify remaining layers
sst_dates <- time(sst_clipped)
print(sst_dates)  # Check if dates match

# Average SST over Niño 3.4 region
sst_ts <- global(sst_clipped, fun = "mean", na.rm = TRUE)
sst_ts <- data.frame(
  date = time(sst_clipped),
  sst_anomaly = sst_ts$mean
)

# Average rainfall over Khartoum
rain_ts <- global(rain_clipped, fun = "mean", na.rm = TRUE)
rain_ts <- data.frame(
  date = time(rain_clipped),
  rainfall = rain_ts$mean
)

# Merge by date
merged_data <- merge(sst_ts, rain_ts, by = "date")

correlation <- cor(merged_data$sst_anomaly, merged_data$rainfall, use = "complete.obs")
print(paste("Correlation Coefficient:", round(correlation, 2)))

# Extract year from date
sst_ts$year <- format(sst_ts$date, "%Y")

# Aggregate monthly SST anomalies to yearly means
sst_yearly <- aggregate(sst_anomaly ~ year, data = sst_ts, FUN = mean, na.rm = TRUE)

# Extract year from date
rain_ts$year <- format(rain_ts$date, "%Y")

# Aggregate monthly rainfall to yearly totals (or means)
rain_yearly <- aggregate(rainfall ~ year, data = rain_ts, FUN = sum, na.rm = TRUE)  # Use sum for total annual rainfall
# OR
rain_yearly <- aggregate(rainfall ~ year, data = rain_ts, FUN = mean, na.rm = TRUE)  # Use mean for annual average

merged_yearly <- merge(sst_yearly, rain_yearly, by = "year")
colnames(merged_yearly) <- c("year", "sst_anomaly", "rainfall")  # Rename columns

correlation_yearly <- cor(merged_yearly$sst_anomaly, merged_yearly$rainfall, use = "complete.obs")
print(paste("Yearly Correlation Coefficient:", round(correlation_yearly, 2)))

ggplot(merged_yearly, aes(x = sst_anomaly, y = rainfall)) +
  geom_point(size = 3, color = "darkblue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(
    title = paste("Yearly SST Anomaly vs. Rainfall in Khartoum (r =", round(correlation_yearly, 2), ")"),
    x = "SST Anomaly (°C)",
    y = "Rainfall (mm)"
  ) +
  theme_minimal()



ggplot(merged_yearly, aes(x = as.numeric(year))) +
  geom_line(aes(y = sst_anomaly, color = "SST Anomaly"), linewidth = 1) +
  geom_line(aes(y = rainfall / 10, color = "Rainfall (scaled)"), linewidth = 1) +  # Scale rainfall for visibility
  scale_y_continuous(
    name = "SST Anomaly (°C)",
    sec.axis = sec_axis(~ . * 10, name = "Rainfall (mm)")
  ) +
  scale_color_manual(values = c("SST Anomaly" = "red", "Rainfall (scaled)" = "blue")) +
  labs(
    title = "Yearly Trends: SST Anomaly(Nino3.4) vs. Rainfall (Khartoum)",
    x = "Year"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())


  #### create animated plots ____-----------------for the SST in Nino3.4-------------------

library(dplyr)
library(ggplot2)
library(gganimate)
library(gifski)

# 1. Aggregate to yearly SST anomalies (already done in your code)
# sst_yearly <- aggregate(sst_anomaly ~ year, data = sst_ts, FUN = mean, na.rm = TRUE)

# 2. Prepare yearly data (ensure chronological order)
merged_yearly <- merged_yearly %>%
  mutate(year_num = as.numeric(year)) %>%  # Convert year to numeric for plotting
  arrange(year_num)  # Ensure ascending order

# 3. Create animated line plot
animated_yearly <- ggplot(merged_yearly, aes(x = year_num, y = sst_anomaly)) +
  geom_line(color = "#e52620", linewidth = 1) +
  geom_point(color = "#f59d22", size = 3) +
  labs(
    title = "Yearly SST Anomaly in Niño 3.4 Region",
    x = "Year",
    y = "SST Anomaly (°C)",
    caption = "Data: NOAA OI SST v2.1"
  ) +
  scale_x_continuous(breaks = seq(min(merged_yearly$year_num), max(merged_yearly$year_num), by = 5)) +  # Ticks every 5 years
  theme_minimal() +
  transition_reveal(year_num) +  # Animate year-by-year
  shadow_mark(color = "gray70", size = 1)  # Keep past points visible

# 4. Render animation
gif <- gganimate::animate(
  animated_yearly,
  nframes = nrow(merged_yearly) * 2,  # 2 frames per year
  fps = 5,
  width = 1000,
  height = 600,
  renderer = gifski_renderer(),
  end_pause = 15
)

# 5. Save
gganimate::anim_save("yearly_sst_animation.gif", animation = gif)

  #### create animated plots ____-----------------for the Rainfall-------------------
library(dplyr)
library(ggplot2)
library(gganimate)
library(gifski)

# 1. Prepare yearly rainfall data
merged_yearly <- merged_yearly %>%
  mutate(year_num = as.numeric(year)) %>%
  arrange(year_num)  # Ensure chronological order

# 2. Create animated rainfall plot
animated_rainfall <- ggplot(merged_yearly, aes(x = year_num, y = rainfall)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "darkorange", size = 3) +
  labs(
    title = "Yearly Rainfall in Khartoum",
    subtitle = "Year: {frame_along}",
    x = "Year",
    y = "Rainfall (mm)",
    caption = "Data: Sudan Meteorological Authority"
  ) +
  scale_x_continuous(breaks = seq(min(merged_yearly$year_num), max(merged_yearly$year_num), by = 5)) +
  theme_minimal() +
  transition_reveal(year_num) +
  shadow_mark(color = "gray70", size = 1)  # Show historical points

# 3. Render animation
rain_gif <- gganimate::animate(
  animated_rainfall,
  nframes = nrow(merged_yearly) * 2,  # 2 frames per year
  fps = 5,
  width = 1000,
  height = 600,
  renderer = gifski_renderer(),
  end_pause = 15
)

# 4. Save
gganimate::anim_save("yearly_rainfall_animation.gif", animation = rain_gif)



  #### create animated plots ____-----------------for the conaind charts-------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(gganimate)
library(gifski)

# 1. Reshape data into long format
merged_yearly_long <- merged_yearly %>%
  pivot_longer(
    cols = c(sst_anomaly, rainfall),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = factor(
      variable,
      levels = c("sst_anomaly", "rainfall"),
      labels = c("SST Anomaly (°C)", "Rainfall (mm)")
    )
  )

# 2. Create combined plot with custom colors
combined_animation <- ggplot(merged_yearly_long, aes(x = as.numeric(year), y = value, color = variable)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  facet_wrap(
    ~variable,
    nrow = 2,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = c(
      "SST Anomaly (°C)" = "#d83417",  # Steel blue
      "Rainfall (mm)" = "#1244ba"       # Orange
    )
  ) +
  scale_y_continuous(expand = c(0.1, 0.1)) +  # Add padding to y-axis
  labs(
    title = "Yearly Climate Trends: Niño 3.4 vs. Khartoum",
    x = "Year",
    caption = "By : Mo Anwar"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    strip.text = element_text(size = 14, face = "bold", color = "white"),  # Facet titles
    strip.background = element_rect(fill = "#333333"),                     # Facet title background
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"  # Remove legend (facets already labeled)
  ) +
  transition_reveal(as.numeric(year)) +
  shadow_mark(color = "gray60", size = 1, alpha = 0.5)

# 3. Render animation
gif <- gganimate::animate(
  combined_animation,
  nframes = 150,
  fps = 10,
  width = 1000,
  height = 800,
  renderer = gifski_renderer(),
  end_pause = 20
)

# 4. Save
gganimate::anim_save("styled_combined_animation.gif", animation = gif)


#___________________________________ animating SST map
# Load required libraries
# Install missing packages first (if needed)


# Load libraries
library(terra)
library(tidyverse)
library(gganimate)
library(gifski)

# 1. Load and preprocess SST data
sst_file <- "D:/Geo x/SST and rainfall/sst.mon.anom.nc"
sst_rast <- rast(sst_file) %>% rotate()

# 2. Aggregate to yearly averages
time_values <- time(sst_rast)
years <- format(time_values, "%Y")
sst_yearly <- tapp(sst_rast, index = years, fun = mean, na.rm = TRUE)

# 3. Prepare dataframe
sst_df <- as.data.frame(sst_yearly, xy = TRUE, na.rm = TRUE) %>%
  pivot_longer(
    cols = -c(x, y),
    names_to = "Year",
    values_to = "SST_Anomaly"
  ) %>%
  mutate(Year = as.numeric(gsub("^X", "", Year))) %>%
  filter(Year >= 1981 & Year <= 2021)

# 4. Define Niño 3.4 region coordinates
nino_box <- data.frame(
  xmin = -170, xmax = -120,
  ymin = -5, ymax = 5
)

# 5. Create SST Map Animation
sst_animation <- ggplot(sst_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = SST_Anomaly)) +
  geom_rect(
    data = nino_box,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    color = "#000000", fill = NA, linewidth = 1.2, linetype = "dashed",
    inherit.aes = FALSE
  ) +
  annotate("text", x = -145, y = 7, label = "Niño 3.4 Region",
           color = "#000000", size = 4, fontface = "bold") +
  borders("world", colour = "black", size = 0.3) +
  scale_fill_gradientn(
    colors = c("#00441b", "#238b45", "#74c477", "#ffff94",
               "#fd8c3c", "#e31a1c", "#800026"),
    values = scales::rescale(c(-3, -1.5, -0.5, 0, 0.5, 1.5, 3)),
    limits = c(-3, 3),
    name = "SST Anomaly (°C)"
  ) +
  coord_fixed(ratio = 1.3) +
  labs(title = "Global Sea Surface Temperature Anomalies") +
  theme_void() +
  transition_time(Year) +  # Critical for animation
  ease_aes("linear")

# 6. Render and save animation
sst_gif <- gganimate::animate(
  sst_animation,
  nframes = length(unique(sst_df$Year)) * 3,
  fps = 8,
  width = 1500,
  height = 900,
  renderer = gifski_renderer(),  # Fixed syntax
  end_pause = 25
)

gganimate::anim_save("sst_map_animation.gif", sst_gif)

  #### create animated plots ____-----------------for SST map and SST chartline -------------------

library(magick)

# Define desired number of frames and fps
desired_frames <- 123  # e.g., 41 years * 3 frames per year
desired_fps <- 5       # slower pace

# Read in the pre-rendered animations and reanimate them to have the same frame count and fps
chart_gif <- image_read("yearly_sst_animation.gif")
chart_gif <- image_animate(chart_gif, nframes = desired_frames, fps = desired_fps)

map_gif <- image_read("sst_map_animation.gif")
map_gif <- image_animate(map_gif, nframes = desired_frames, fps = desired_fps)

# Ensure both GIFs now have the same number of frames
nframes <- min(length(chart_gif), length(map_gif))
chart_gif <- chart_gif[1:nframes]
map_gif <- map_gif[1:nframes]

# Define the full range of years (1981-2021) and calculate frames per year
years <- 1981:2021
n_years <- length(years)
frames_per_year <- nframes / n_years  # frames corresponding to one year

# Prepare a list to store combined frames
combined_frames <- vector("list", nframes)

for(i in 1:nframes) {
  
  # Adjust chart frame: add a vertical margin (30px) to push it down
  chart_frame <- chart_gif[i]
  margin <- image_blank(
    width = image_info(chart_frame)$width, 
    height = 30, 
    color = "white"
  )
  chart_frame_down <- image_append(c(margin, chart_frame), stack = TRUE)
  
  # Resize the chart frame (with margin) to have the same height as the map frame
  map_height <- image_info(map_gif[i])$height
  chart_frame_down <- image_resize(chart_frame_down, paste0("x", map_height))
  
  # Combine the adjusted chart (left) with the map (right) side-by-side
  combined_side <- image_append(c(chart_frame_down, map_gif[i]), stack = FALSE)
  
  # Compute the current year for this frame (assumes even distribution)
  current_year <- 1980 + ceiling(i / frames_per_year)
  
  # Create a header (blank image) with the same width as the combined frame and 100px height
  header <- image_blank(
    width = image_info(combined_side)$width, 
    height = 100, 
    color = "white"
  )
  
  # Annotate the header with a common title and dynamic year label
  header <- image_annotate(
    header, 
    text = paste("Global Sea Surface Temperature Anomalies  |  Year:", current_year),
    size = 40,
    font = "sans",
    color = "black",
    gravity = "center"
  )
  
  # Stack the header above the combined frame
  final_frame <- image_append(c(header, combined_side), stack = TRUE)
  
  # Save this final frame
  combined_frames[[i]] <- final_frame
}

# Join all frames into a single animated image
final_animation <- image_join(combined_frames)

# Animate the final image at the desired fps (which is a factor of 100)
final_animation <- image_animate(final_animation, fps = desired_fps)

# Save the final combined animation
image_write(final_animation, path = "combined_animation2.gif")
