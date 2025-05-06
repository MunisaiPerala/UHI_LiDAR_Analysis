# UHI_LiDAR_Analysis
Repository for Urban Heat Island (UHI) Analysis using LiDAR and R.

# ️ Urban Heat Island (UHI) Analysis Using LiDAR

##  Project Overview
This project explores how vegetation structure (tree canopy cover) and terrain (elevation) affect surface temperature patterns in urban environments. The analysis uses LiDAR data, Landsat 8 thermal imagery, and NLCD land cover.

## Tools Used
- R programming language
- Packages: `terra`, `sf`, `lidR`, `ggplot2`, `scatterplot3d`
- Software: RStudio

##  Folder Structure
GIT H/37/ 
├── data_raw/       # Input rasters and .las files
├── outputs/        # Plots, processed rasters, CSVs
├── scripts/        # R scripts for processing and modeling
├── report/         # Final project report DOCX
├── README.md       # This file

## ️ Workflow Summary
- Harmonize CRS and extent of temperature and canopy datasets
- Create 1 km² spatial grids over the study area
- Extract mean temperature, tree canopy %, and elevation per grid cell
- Run linear and multiple regression models to explore UHI patterns
- Visualize results using 2D, Histograms and 3D plots

##  Key Findings
- Tree cover and elevation did not significantly explain variation in surface temperature
- However, visualizations showed spatial patterns consistent with the UHI effect
- The workflow demonstrates integration of multiple spatial datasets in R

## Author
**Muni Sai Ganesha Perala**  
Master’s in Bioinformatics, Saint Louis University
