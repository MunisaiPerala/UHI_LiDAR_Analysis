# ---------------------------------------------
# Urban Heat Island (UHI) Analysis Project
# Author: Muni Sai Ganesha Perala
# Description: Analyzing UHI patterns using Landsat, NLCD, and LiDAR in R
# Date: May 2025
# ---------------------------------------------

# ----------------------------
# 1. Load Required Packages
# ----------------------------
install.packages(c("terra", "sf", "lidR", "ggplot2", "scatterplot3d"))
library(terra)
library(sf)
library(lidR)
library(ggplot2)
library(scatterplot3d)

# ----------------------------
# 2. Load and Preprocess Data
# ----------------------------

# Load Landsat 8 Band 10 - Thermal Infrared (surface temperature)
temp_raster <- rast("data_raw/LC09_L2SP_B10.TIF")

# Load 2021 Tree Canopy Cover raster from NLCD
tree_cover <- rast("data_raw/nlcd_tcc_conus_2021.tiff")

# Reproject tree canopy raster to match temperature CRS
tree_aligned <- project(tree_cover, temp_raster)

# Crop both rasters to shared extent
common_extent <- intersect(ext(temp_raster), ext(tree_aligned))
temp_crop <- crop(temp_raster, common_extent)
tree_crop <- crop(tree_aligned, common_extent)

# Quick check plots
plot(temp_crop, main = "Cropped Surface Temp (°C)")
plot(tree_crop, main = "Cropped Tree Cover (%)")

# ----------------------------
# 3. Create Grid and Extract Zonal Statistics
# ----------------------------

# Generate 1 km grid
grid_size <- 1000  # meters
study_ext <- ext(temp_crop)
grid_vect <- rast(extent = study_ext, resolution = grid_size)
grid_vect_poly <- as.polygons(grid_vect)
grid_vect_poly$ID <- 1:nrow(grid_vect_poly)

# Extract average values per grid cell
mean_temp <- extract(temp_crop, grid_vect_poly, fun = mean, na.rm = TRUE)
mean_tree <- extract(tree_crop, grid_vect_poly, fun = mean, na.rm = TRUE)

# Merge into dataframe
grid_data <- data.frame(
  ID = mean_temp$ID,
  mean_temp_c = mean_temp[, 2],
  mean_tree_cover = mean_tree[, 2]
)
grid_data_clean <- na.omit(grid_data)

# ----------------------------
# 4. LiDAR-Based Terrain and Canopy Modeling
# ----------------------------

# Load LiDAR point cloud (.las file)
las <- readLAS("data_raw/3809022_SE.las")

# Filter for ground points and create Digital Terrain Model (DTM)
las_ground <- filter_poi(las, Classification == 2)
dtm <- rasterize_terrain(las_ground, res = 10, algorithm = knnidw(k = 10))
plot(dtm, main = "Digital Terrain Model (DTM)")
plot(grid_vect_poly, add = TRUE, border = "red")

# Align and crop DTM
dtm_aligned <- project(dtm, temp_crop)
dtm_crop <- crop(dtm_aligned, ext(grid_vect_poly))

# Extract elevation per grid cell
mean_elev <- extract(dtm_crop, grid_vect_poly, fun = mean, na.rm = TRUE)
elev_df_clean <- na.omit(data.frame(ID = mean_elev$ID, mean_elevation = mean_elev[, 2]))

# Merge with main grid dataset
grid_data_final <- merge(grid_data_clean, elev_df_clean, by = "ID")

# ----------------------------
# 5. Statistical Modeling
# ----------------------------


# Model 1: Tree Canopy Cover Only

# The first linear model tested whether tree canopy cover alone could explain 
# variation in surface temperature. The results showed a very low R² (~0.003) 
# and a high p-value (0.7157), indicating no significant relationship. 
# Interestingly, the coefficient was positive (+17.55), which was opposite 
# of what was expected — but it's not meaningful due to the lack of significance.

uhi_model <- lm(mean_temp_c ~ mean_tree_cover, data = grid_data_clean)
summary(uhi_model)

# Model 2: Temp ~ Tree Canopy + Elevation

# The second model included elevation as an additional predictor.
# The R² increased slightly to ~0.034, but neither tree canopy (p = 0.617) 
# nor elevation (p = 0.447) were statistically significant.
# Tree cover had a negative coefficient (as expected), but the model still 
# explains very little variation in temperature.

uhi_model2 <- lm(mean_temp_c ~ mean_tree_cover + mean_elevation, data = grid_data_final)
summary(uhi_model2)

# Save final dataset and model output
write.csv(grid_data_final, "outputs/uhi_grid_data_final.csv", row.names = FALSE)
sink("outputs/uhi_model_summary.txt")
summary(uhi_model2)
sink()

# Overall, this analysis suggests that tree cover and elevation alone may not be 
# strong predictors of surface temperature in this dataset. Other factors such as 
# impervious surfaces, time of image capture, or urban material type may be more 
# influential in driving Urban Heat Island effects.

# ----------------------------
# 6. Visualizations
# ----------------------------

# 2D Scatterplot with regression line
ggplot(grid_data_clean, aes(x = mean_tree_cover, y = mean_temp_c)) +
  geom_point(color = "darkgreen", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(title = "Tree Cover vs Surface Temp", x = "Tree Canopy (%)", y = "Temp (°C)") +
  theme_minimal()

# 3D Scatterplot of Tree Cover, Elevation, and Temperature
scatterplot3d(
  x = grid_data_final$mean_tree_cover,
  y = grid_data_final$mean_elevation,
  z = grid_data_final$mean_temp_c,
  color = "darkblue",
  pch = 19,
  xlab = "Tree Canopy (%)",
  ylab = "Elevation (m)",
  zlab = "Temperature (°C)",
  main = "3D View: Tree Cover, Elevation, Temp"
)

# ----------------------------
# 7. Canopy Height Modeling (Optional)
# ----------------------------

# Digital Surface Model from all points
dsm <- rasterize_canopy(las, res = 10, algorithm = p2r())
chm <- dsm - dtm
plot(chm, main = "Canopy Height Model (CHM)", col = terrain.colors(20))

# CHM Histogram
hist(values(chm)[!is.na(values(chm))], breaks = 20, col = "forestgreen",
     main = "Distribution of Canopy Heights", xlab = "Height (m)")

# Save raster
writeRaster(chm, "outputs/chm.tif", overwrite = TRUE)