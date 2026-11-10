
# ==============================================================================
# SUMMARISE THE MONTHLY CLIMATOLOGIES OF CHELSA TO COMPARE WITH ANTTERRCLIM
# ==============================================================================

# Takes every monthly CHELSA .tif file which has been downloaded in script 2a and summarises them into monthly climatologies following the methods of script 1b

# It also calculates the mean annual temperature, mean austral summer temperature and mean austral winter temperature

# This was run on a local computer 

# The outputs are saved in the "Data/CHELSA" folder

# Single job and the rasters are loaded per (month, year) 

# No reprojection and no ice-free masking are applied in this script.


library(dplyr)
library(purrr)
library(terra)
library(here)


# Set inpath -------------------------------------------------------------

inpath <- here("Data/CHELSA")

# Set outpath -------------------------------------------------------------

outpath <- here("Data/CHELSA/comparison")


###################################
# HISTORICAL (1981 - 2010) ------------------------------------------------
###################################

hist_files <- list.files(inpath, pattern = "CHELSA_tas_\\d{2}_1981-2010_V\\.2\\.1\\.tif$", 
                         full.names = TRUE) %>% 
  sort() # ensures months 01-12 order

chelsa_hist <- rast(hist_files)

# Name layers by month
names(chelsa_hist) <- sprintf("tas_%02d", 1:12)


# Annual temp -------------------------------------------------------------

mean_annual_temp_K <- app(chelsa_hist, mean, na.rm = TRUE)
# Convert from Kelvin to Celsius
mean_annual_temp <- mean_annual_temp_K - 273.15

writeRaster(mean_annual_temp, file.path(outpath, "Mean_Annual_Temperature_HISTORICAL_1981_2010.tif"), overwrite = T)



##########################
# Summer temp (DJF) --------------------------------------------------------
##########################

chelsa_hist_DJF <- chelsa_hist[[c(1, 2, 12)]]

mean_summer_temp_K <- app(chelsa_hist_DJF, mean, na.rm = TRUE)
# Convert from Kelvin to Celsius
mean_summer_temp <- mean_summer_temp_K - 273.15

writeRaster(mean_summer_temp, file.path(outpath, "Mean_Summer_Temperature_HISTORICAL_1981_2010.tif"), overwrite = T)


######################
# Winter temp (DJF) --------------------------------------------------------
######################

chelsa_hist_JJA <- chelsa_hist[[c(6, 7, 8)]]

mean_winter_temp_K <- app(chelsa_hist_JJA, mean, na.rm = TRUE)
# Convert from Kelvin to Celsius
mean_winter_temp <- mean_winter_temp_K - 273.15

writeRaster(mean_winter_temp, file.path(outpath, "Mean_Winter_Temperature_HISTORICAL_1981_2010.tif"), overwrite = T)


###################################
# FUTURE (2071 - 2100) ------------------------------------------------
###################################

future_models  <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0",  "ukesm1-0-ll")

model <- future_models[job_index]

for(model in future_models){
  
  print(paste0("Running model: ", model))
  
  future_files <- list.files(inpath, pattern = paste0("CHELSA_", model, "_r1i1p1f1_w5e5_ssp370_tas_\\d{2}_2071-2100_V\\.2\\.1\\.tif$"),
                             full.names = TRUE)
  
  chelsa_future <- rast(future_files)
  
  # Name layers by month
  names(chelsa_future) <- sprintf("tas_%02d", 1:12)
  
  # Annual temp -------------------------------------------------------------
  mean_annual_temp_K <- app(chelsa_future, mean, na.rm = TRUE)
  # Convert from Kelvin to Celsius
  mean_annual_temp <- mean_annual_temp_K - 273.15
  writeRaster(mean_annual_temp, file.path(outpath, paste0("Mean_Annual_Temperature_FUTURE_", model, "_2071_2100.tif")), overwrite = T)
  
  
  ##########################
  # Summer temp (DJF) --------------------------------------------------------
  ##########################
  chelsa_future_DJF <- chelsa_future[[c(1, 2, 12)]]
  mean_summer_temp_K <- app(chelsa_future_DJF, mean, na.rm = TRUE)
  # Convert from Kelvin to Celsius
  mean_summer_temp <- mean_summer_temp_K - 273.15
  writeRaster(mean_summer_temp, file.path(outpath, paste0("Mean_Summer_Temperature_FUTURE_", model, "_2071_2100.tif")), overwrite = T)

  
  ######################
  # Winter temp (JJA) --------------------------------------------------------
  ######################
  chelsa_future_JJA <- chelsa_future[[c(6, 7, 8)]]
  mean_winter_temp_K <- app(chelsa_future_JJA, mean, na.rm = TRUE)
  # Convert from Kelvin to Celsius
  mean_winter_temp <- mean_winter_temp_K - 273.15
  writeRaster(mean_winter_temp, file.path(outpath, paste0("Mean_Winter_Temperature_FUTURE_", model, "_2071_2100.tif")), overwrite = T)
  
  
}

