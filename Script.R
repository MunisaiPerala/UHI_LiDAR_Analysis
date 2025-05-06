# In this project, I chose to focus on three main spatial variables that are commonly associated with the 
# Urban Heat Island (UHI) effect: surface temperature, tree canopy cover, and elevation. My decision was based 
# on both prior research and the types of data I was able to access and work with effectively.

#	Surface temperature was derived from Landsat 8 Band 10 (thermal infrared), a reliable satellite source frequently used for UHI analysis. This allowed me to detect spatial variations in heat across my study area.
#	Tree canopy cover data came from the 2021 NLCD Tree Canopy layer, which provides percent cover information. I included this because vegetation is known to provide cooling through shade and evapotranspiration, and I wanted to test that hypothesis in an urban context.
# Elevation was extracted from high-resolution LiDAR point cloud data, used to build a Digital Terrain Model (DTM). I added this variable because land elevation can influence air flow, drainage, and heat retention—factors that might affect urban temperature patterns.

# For the analysis, I used linear regression because it is simple, interpretable, and suited to exploring 
# relationships between continuous variables. This approach allowed me to evaluate how well tree cover and 
# elevation explained variations in surface temperature across a spatial grid.

# Overall, my analytical choices were grounded in environmental theory, data availability, and the desire to keep the methods approachable and reproducible, especially as I’m still developing my R and GIS skills.



install.packages("terra")
library(terra)
# Load temp data
temp_raster <- rast("data_raw/LC09_L2SP_B10.TIF")

# Load Tree Canopy Cover raster
tree_cover <- rast("data_raw/nlcd_tcc_conus_2021.tiff")

# Quick visualization
plot(tree_cover, main = "Tree Canopy Cover (%)")


# Harmonize CRS and Extent
# Reproject tree_cover to match temperature raster
tree_aligned <- project(tree_cover, temp_raster)

# Crop to shared extent
common_extent <- intersect(ext(temp_raster), ext(tree_aligned))
temp_crop <- crop(temp_raster, common_extent)
tree_crop <- crop(tree_aligned, common_extent)

# Quick check
plot(temp_crop, main = "Cropped Surface Temp (°C)")
plot(tree_crop, main = "Cropped Tree Cover (%)")



#Installing Package 
install.packages("sf")

library(terra)  
library(sf)     # for spatial vector operations


-------------------------------------------------
  # Step 1: Set grid size (1 km = 1000 meters)
  grid_size <- 1000  

# Step 2: Get extent of cropped raster
study_ext <- ext(temp_crop)

# Step 3: Create a regular grid over the extent
grid_vect <- rast(extent = study_ext, resolution = grid_size)

# Step 4: Convert to polygons
grid_vect_poly <- as.polygons(grid_vect)

# Step 5: Optional - assign IDs to each grid cell
grid_vect_poly$ID <- 1:nrow(grid_vect_poly)

# Step 6: Visualize
plot(grid_vect_poly, main = "1 km Grid over Study Area")


----------------------------------------------------------
  # Extract mean surface temperature (°C) per grid cell
  mean_temp <- extract(temp_crop, grid_vect_poly, fun = mean, na.rm = TRUE)

# Extract mean tree canopy cover (%) per grid cell
mean_tree <- extract(tree_crop, grid_vect_poly, fun = mean, na.rm = TRUE)


# Combine into a clean data frame
grid_data <- data.frame(
  ID = mean_temp$ID,
  mean_temp_c = mean_temp[, 2],
  mean_tree_cover = mean_tree[, 2]
)

head(grid_data)

# Remove rows with missing values
grid_data_clean <- na.omit(grid_data)

# Check cleaned result
head(grid_data_clean)

install.packages("ggplot2")
library(ggplot2)

ggplot(grid_data_clean, aes(x = mean_tree_cover, y = mean_temp_c)) +
  geom_point(color = "darkgreen", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    title = "UHI Effect: Tree Cover vs Surface Temperature",
    x = "Mean Tree Canopy Cover (%)",
    y = "Mean Surface Temperature (°C)"
  ) 
  theme_minimal()
  
  
  # lm() function to see how tree canopy cover affects temperature.
  # Fit a linear model: Temperature ~ Tree Cover
  uhi_model <- lm(mean_temp_c ~ mean_tree_cover, data = grid_data_clean)
  
  # View model summary
  summary(uhi_model)
  
  
----------------------------------------------  
  # LIDAR
  
  # Install and Loading {lidR} Package
  install.packages("lidR")   
  
  #Load the LiDAR Point Cloud
  # Load the .las file
  las <- readLAS("data_raw/3809022_SE.las")
  
  # Check if it loaded
  summary(las)
  library(lidR)              
  
  #Build a Digital Terrain Model (DTM)
  # Filter only ground points (Classification == 2)
  las_ground <- filter_poi(las, Classification == 2)
  
  # Create a Digital Terrain Model (10m resolution)
  dtm <- rasterize_terrain(las_ground, res = 10, algorithm = knnidw(k = 10))
  
  # Plot the result
  plot(dtm, main = "Digital Terrain Model (DTM)")
  

  # Plot DTM
  plot(dtm, main = "DTM")
  
    # Overlay the grid
  plot(grid_vect_poly, add = TRUE, border = "red")
  
 
  
  
   # Align DTM with temp_crop raster (so grid and DTM match)
  dtm_aligned <- project(dtm, temp_crop)
  
  # Crop DTM to match extent of your grid
  dtm_crop <- crop(dtm_aligned, ext(grid_vect_poly))
  
  # Extract mean elevation per grid cell
  mean_elev <- extract(dtm_crop, grid_vect_poly, fun = mean, na.rm = TRUE)
  
  # Create elevation data frame
  elev_df <- data.frame(
    ID = mean_elev$ID,
    mean_elevation = mean_elev[, 2]
  )
  
  # Clean out any NA values
  elev_df_clean <- na.omit(elev_df)
  
  # Prepare base data (only non-NA grid cells from earlier)
  grid_base <- grid_data[, c("ID", "mean_temp_c", "mean_tree_cover")]
  grid_base <- na.omit(grid_base)
  
  # Merge elevation with existing tree + temperature values
  grid_data_final <- merge(grid_base, elev_df_clean, by = "ID")
  
  # Check the final table
  head(grid_data_final)
 
  # Final regression model including elevation
  uhi_model2 <- lm(mean_temp_c ~ mean_tree_cover + mean_elevation, data = grid_data_final)
  
  # View the summary of the model
  summary(uhi_model2)
  
  
  write.csv(grid_data_final, "outputs/uhi_grid_data_final.csv", row.names = FALSE)
  
 ------------------------------------------------------------------------------- 
  # Canopy Height Model
  # Step 1: Create DSM (Digital Surface Model — includes treetops/buildings)
  dsm <- rasterize_canopy(las, res = 10, algorithm = p2r())
  
  # Step 2: Calculate CHM (Canopy Height Model = DSM - DTM)
  chm <- dsm - dtm
  
  # Step 3: Plot CHM
  plot(chm, main = "Canopy Height Model (CHM)", col = terrain.colors(20))
  
  
  # Convert raster to values for plotting
  chm_values <- values(chm)
  chm_values <- chm_values[!is.na(chm_values)]  # remove NAs
  
  hist(chm_values,
       breaks = 20,
       col = "forestgreen",
       main = "Distribution of Canopy Heights",
       xlab = "Height in meters")
  
  # Install package 
  install.packages("scatterplot3d")
  
  # Load the package
  library(scatterplot3d)
  
  # 3D Scatterplot of Tree Cover, Elevation, and Temperature
  scatterplot3d(
    x = grid_data_final$mean_tree_cover,
    y = grid_data_final$mean_elevation,
    z = grid_data_final$mean_temp_c,
    color = "darkblue",
    pch = 19,
    xlab = "Tree Canopy Cover (%)",
    ylab = "Elevation (m)",
    zlab = "Surface Temperature (°C)",
    main = "3D View of Urban Heat Island Influences"
  )
  
  # Interpretation of Regression Results
  
  # Model 1: Tree Canopy Cover Only
  # The first linear model tested whether tree canopy cover alone could explain 
  # variation in surface temperature. The results showed a very low R² (~0.003) 
  # and a high p-value (0.7157), indicating no significant relationship. 
  # Interestingly, the coefficient was positive (+17.55), which was opposite 
  # of what was expected — but it's not meaningful due to the lack of significance.
  
  uhi_model <- lm(mean_temp_c ~ mean_tree_cover, data = grid_data_clean)
  
  # Model 2: Tree Canopy + Elevation
  
  # The second model included elevation as an additional predictor.
  # The R² increased slightly to ~0.034, but neither tree canopy (p = 0.617) 
  # nor elevation (p = 0.447) were statistically significant.
  # Tree cover had a negative coefficient (as expected), but the model still 
  # explains very little variation in temperature.
  
  uhi_model2 <- lm(mean_temp_c ~ mean_tree_cover + mean_elevation, data = grid_data_final)
  
  
  # Overall, this analysis suggests that tree cover and elevation alone may not be 
  # strong predictors of surface temperature in this dataset. Other factors such as 
  # impervious surfaces, time of image capture, or urban material type may be more 
  # influential in driving Urban Heat Island effects.
  # -------------------------------------------------------------------
  
  
  # Save CHM raster
  writeRaster(chm, "outputs/chm.tif", overwrite = TRUE)
  
  # Save model summary
  sink("outputs/uhi_model_summary.txt")
  summary(uhi_model2)
  sink()
  # ----------------------------------------------------------------
  
 