
# ==============================================================================
# SUMMARISE THE RAW CLIMATE VARIABLES INTO ANNUAL BIOCLIMATIC VARIABLES
# ==============================================================================

# This script loads the climate variables summarise with scripts 1b. in a monthly climatology format and summarises them to bioclimatic variables (12 months each of TasMin, TasMax, and Total Precipitation, per model and period) and runs them through predicts::bcvars() to get the standard 19 WorldClim-style bioclimatic variables, saved as one .tif per BIO layer (BIO1, BIO2, ... BIO19) per model/period.

# DEPENDENCY: this needs scripts 1b. "temp_minmax" variable to have been run (it produces the TasMin/TasMax monthly climatologies) and "total_precip" to have been run (monthly precipitation totals). If either is missing for a given model/period, this stops with a message telling you exactly which file(s) are missing.

# Designed to be run on a High-Performance Computer with an array of jobs but could probably be run on a local computer

# One job = one model (all periods for that model are done in the same job).

# The script creates a folder directory structured as: here("Data/Environmental_predictors/PolarRes26", model) where the monthly climatologies are saved as .tif files

# The outputs are still in their native extent and resolution. The next step will harmonise them to a common grid.

# * NOTE - each job processes hundreds of rasters and took from 12-40 hours on a high-performance computer

# * Also NOTE - a couple of these variables will fail for HCLIM_CESM2,  HCLIM_MPI_ESM1 and HCLIM_ERA5 (ones involving siconca, rsus, snm). That's because they were downloaded separately (in the _Additional scripts) and will be summarised in separate scripts due to slightly different data structures (1b_Summarising_Climate_Vars_HCLIM_CESM2_Additional.R etc.)


# ---- 0. Setup ------------------------------------------------------------------

lib_loc <- paste(getwd(), "/r_lib_new", sep = "")

library(terra)
library(here)
library(predicts, lib.loc = lib_loc)

# ---- 1. Configuration -----------------------------------------------------------

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5", "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5", "MetUM_ERA5")

years_hist   <- seq(1995, 2014, by = 1)
years_mid    <- seq(2041, 2060, by = 1)
years_future <- seq(2081, 2100, by = 1)

input_base  <- here("Data/Environmental_predictors/PolarRes26")
output_base <- here("Data/Environmental_predictors/PolarRes26_Bioclim")

# One job = one model. job_index (1-7) picks a row from `models`.
args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])


model <- models[job_index]
message("Job ", job_index, "/", length(models), " -> Model: ", model)

# Just the period labels + filename range strings needed to find the
# already-computed monthly climatology outputs -- doesn't need
# years/scenario/RACMO-clipping info the way the main summarising script
# does, since no daily data is read here.
#
# ERA5-driven models (HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5) only ever have a
# HISTORICAL set of monthly climatology files

is_era5 <- endsWith(model, "_ERA5")

periods <- if (is_era5) {
  list(HISTORICAL = paste(min(years_hist), max(years_hist), sep = "_"))
} else {
  list(
    HISTORICAL = paste(min(years_hist),   max(years_hist),   sep = "_"),
    MID        = paste(min(years_mid),    max(years_mid),    sep = "_"),
    FUTURE     = paste(min(years_future), max(years_future), sep = "_")
  )
}

model_in_dir <- file.path(input_base, model)
if (!dir.exists(model_in_dir)) {
  stop("Model directory not found: ", model_in_dir)
}

model_out_dir <- file.path(output_base, model)
dir.create(model_out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Helper: load one variable's 12 monthly climatology files, in order -----

load_monthly_climatology <- function(model_dir, name_pattern, period_name, range_label) {
  filenames <- sprintf(name_pattern, month.name, period_name, range_label)
  paths <- file.path(model_dir, filenames)
  
  missing <- !file.exists(paths)
  if (any(missing)) {
    stop("Missing monthly climatology file(s) -- has the matching variable ",
         "been run for this model/period yet?\n",
         paste(paths[missing], collapse = "\n"))
  }
  
  # Loaded in month.name (Jan-Dec) order, matching the order `filenames`
  # was built in -- this is what bcvars() requires.
  rast(paths)
}

# ---- 3. Build bioclim variables for this model, every period -------------------

for (period_name in names(periods)) {
  range_label <- periods[[period_name]]
  message("-- bioclim: ", model, " / ", period_name)
  
  tmin <- load_monthly_climatology(model_in_dir, "Climatological_Monthly_Mean_TasMin_%s_%s_%s.tif",
                                   period_name, range_label)
  tmax <- load_monthly_climatology(model_in_dir, "Climatological_Monthly_Mean_TasMax_%s_%s_%s.tif",
                                   period_name, range_label)
  prec <- load_monthly_climatology(model_in_dir, "Climatological_Monthly_Total_Precipitation_%s_%s_%s.tif",
                                   period_name, range_label)
  
  bio <- predicts::bcvars(prec, tmin, tmax, filename = "")
  
  # Save each BIO layer separately (BIO1, BIO2, ... BIO19) rather than one
  # multi-layer file
  for (i in seq_len(nlyr(bio))) {
    out_path <- file.path(model_out_dir, sprintf("BIO%d_%s_%s.tif", i, period_name, range_label))
    writeRaster(bio[[i]], out_path, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  }
  message("  wrote BIO1-BIO", nlyr(bio), " for ", period_name, " to ", model_out_dir)
}

message("Done: ", model)