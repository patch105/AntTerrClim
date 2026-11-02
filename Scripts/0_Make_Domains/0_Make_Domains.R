
###############################################
# Load libraries -----------------------------------------------------------
################################################

library(terra)
library(here)

###############################################
# Study domains -----------------------------------------------------------
################################################

# Set the boundary for our raster grids from Gerrish (2020) 60S boundary:

Ant_extent <- vect(here("Data/add_data_limit_v7.2.shp")) # EPSG:3031
Ant_extent_grid_1km <- rast(Ant_extent, res = 1000, vals = 1)
Ant_extent_grid_10km <- rast(Ant_extent, res = 10000, vals = 1)
Ant_extent_grid_100m <- rast(Ant_extent, res = 100, vals = 1)

writeRaster(Ant_extent_grid_1km, here("Data/Ant_extent_grid_1km.tif"), overwrite =T)
writeRaster(Ant_extent_grid_10km, here("Data/Ant_extent_grid_10km.tif"), overwrite = T)

# Load the two domains for AntTerrEnv variables:

# 1. Antarctic continent, from SCAR Antarctic Digital Database high resolution vector polygons of the Antarctic coastline v7.12 (Gerrish et al., 2026) 
# NOTE: Downloaded the zipped file and unzipped it manually

coast <- vect(here("Data/add_coastline_high_res_polygon_v7_12.shp"))

coast_domain <- mask(Ant_extent_grid_10km, coast, updatevalue = NA, touches = T)

names(coast_domain) <- "coast_domain" 

writeRaster(coast_domain, here("Data/coast_domain.tif"), overwrite = T)


# 1b. Make a ocean mask from the coast - but keeping the maximum extent of ocean in with touches=T
ocean_domain <- mask(Ant_extent_grid_10km, coast, inverse = T, updatevalue = NA, touches = T)

names(ocean_domain) <- "ocean_domain" 

writeRaster(ocean_domain, here("Data/ocean_domain.tif"), overwrite = T)


# 2A. Ice-free areas, from Tóth and Terauds (2023) - polygon version supplied directly by Anikó Tóth

ice_free <- vect(here("Data/rocks_Union_Land.shp"))

ice_free_domain <- mask(Ant_extent_grid_10km, ice_free, updatevalue = NA, touches = T)

names(ice_free_domain) <- "ice_free_domain" 

writeRaster(ice_free_domain, here("Data/ice_free_domain.tif"), overwrite = T)


# Also make a 100m ice-free area version for some later calculations

ice_free_domain_100m <- mask(Ant_extent_grid_100m, ice_free, updatevalue = NA, touches = T)

names(ice_free_domain_100m) <- "ice_free_domain" 

writeRaster(ice_free_domain_100m, here("Data/ice_free_domain_100m.tif"), overwrite = T)


# Also make a 1km ice-free area version for some later calculations

ice_free_domain_1km <- mask(Ant_extent_grid_1km, ice_free, updatevalue = NA, touches = T)

names(ice_free_domain_1km) <- "ice_free_domain" 

writeRaster(ice_free_domain_1km, here("Data/ice_free_domain_1km.tif"), overwrite = T)


# 2B. ** OPTIONAL: FUTURE ice-free areas, from Lee et al. (2017) 
# Downloaded from: http://dx.doi.org/doi:10.4225/15/585216f8703d0 (Lee & Terauds, 2017) 

# Load the future ice-free layer as a shapefile
ice_free_future <- vect(here("Data/AAS_4297_Future_Ice-free_Layers/AAS_4297_Ice_Free_Shapefiles/PS_RCP45_Best_Future_IceFree.shp"))

ice_free_future_domain <- mask(Ant_extent_grid_10km, ice_free_future, updatevalue = NA, touches = T)

names(ice_free_future_domain) <- "ice_free_future_domain" 

writeRaster(ice_free_future_domain, here("Data/ice_free_future_domain.tif"), overwrite = T)


# Also make a 100m future ice-free area version for some later calculations

ice_free_future_domain_100m <- mask(Ant_extent_grid_100m, ice_free_future, updatevalue = NA, touches = T)

names(ice_free_future_domain_100m) <- "ice_free_future_domain" 

writeRaster(ice_free_future_domain_100m, here("Data/ice_free_future_domain_100m.tif"), overwrite = T)



###############################################
# References --------------------------------------------------------------
################################################

# Gerrish, L. (2020). Antarctic Digital Database data limit at 60S (7.2) [Data set]. UK Polar Data Centre, Natural Environment Research Council, UK Research & Innovation. https://doi.org/10.5285/367C992B-55A7-4A90-B972-861E443F95A1. [Accessed 30 June 2026]

# Gerrish, L., Ireland, L., Fretwell, P., Cooper, P., & Skachkova, A. (2026). High resolution vector polygons of the Antarctic coastline 2 (Version 7.12) [Data set]. NERC EDS UK Polar Data Centre. https://doi.org/10.5285/13c4d2f1-8903-4d7f-8977-592121975554 [Accessed 09 June 2026]

# Lee, J.R., Raymond, B., Bracegirdle, T.J., Chadès, I., Fuller, R.A., Shaw, J.D., Terauds, A., 2017. Climate change drives expansion of Antarctic ice-free habitat. Nature 547, 49–54. https://doi.org/10.1038/nature22996

# Lee, J., and Terauds, A. (2017) Projections of Antarctic ice-free areas under two RCP climate change scenarios, Ver. 1, Australian Antarctic Data Centre - doi:10.4225/15/585216f8703d0 [Accessed 10 November 2025]

# Tóth, A. and Terauds, A. (2023). Ice-free Antarctica - A union of rock outcrop layers derived from imagery collected from 1960-2020, Ver. 1, Australian Antarctic Data Centre – https://doi:10.26179/7mnh-j215 [Accessed 04 February 2025]  
