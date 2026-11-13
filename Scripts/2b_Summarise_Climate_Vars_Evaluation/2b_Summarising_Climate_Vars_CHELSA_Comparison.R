# ==============================================================================
# SUMMARISE TEMPERATURE INTO MONTHLY CLIMATOLOGIES FOR MODEL EVALUATION (CHELSA COMPARISON)
# ==============================================================================

# This script loads the daily near-surface air temperature (tas) downloaded with scripts 1a. and summarises it to monthly climatologies, for the evaluation periods 1985-2010 (HISTORICAL) and 2071-2100 (FUTURE).

# It is a trimmed version of the full climatology script (1b.): it computes ONLY temperature (12 monthly climatologies, mean annual, and DJF/JJA seasonal means), across all 7 models. The periods differ from the full script's so that the outputs can be compared against the CHELSA climatologies.

# Designed to be run on a High-Performance Computer with an array of jobs

# Every job is a separate model - specified from the 'models' list below (job index 1-7)

# The script creates a folder directory structured as: here("Data/Environmental_predictors/PolarRes26", model, "comparison") where the monthly climatologies are saved as .tif files. The extra "comparison" subfolder keeps these outputs separate from the main climatologies made by the full script.

# The outputs are still in their native extent and resolution. No resampling or harmonisation is done in this script.

# * NOTE - the historical and future periods are built from different scenario folders depending on the model:
#   - HISTORICAL: scenario "evaluation" for the ERA5-driven models (HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5), scenario "historical" for the rest.
#   - FUTURE: scenario "ssp370", computed only for the non-ERA5 models, since the ERA5-driven models have no ssp series (same logic as the full script).

# * NOTE - each job processes hundreds of rasters and takes from 12-24 hours on a high-performance computer

# * Also NOTE - this script assumes the diagnostics from Script 1 have already been checked (no spillover, no gaps, no duplicate dates, daily data, standard calendar, tas in K).

# ---- 0. Setup ------------------------------------------------------------------

library(terra)
library(here)
library(lubridate)
library(stringr)
library(dplyr)
library(purrr)
library(ncdf4)

# ---- 1. Configuration -----------------------------------------------------------

base_dir <- here("Data/PolarRes26")

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5",
            "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5",
            "MetUM_ERA5")

ssp_scenario <- "ssp370"

years_hist   <- seq(1985, 2010, by = 1)
years_future <- seq(2071, 2100, by = 1)

# One job = one model. args[1] = model index (1-7), matching the order of
# the 'models' list above.
args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])

model <- models[job_index]

# ERA5-driven models (HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5) don't have an ssp
# series at all, and their historical-equivalent files live under an
# "evaluation" folder rather than "historical" -- same logic as the full
# climatology script.
is_era5 <- endsWith(model, "_ERA5")

# RACMO's ssp series currently ends 2099-12-31 (one year short of HCLIM's
# 2100). For RACMO, the years actually used to compute the FUTURE period
# are clipped to 2099, but the output filenames keep the original
# "2071_2100" label (set below) so RACMO and HCLIM outputs stay directly
# comparable.
racmo_future_end_year <- 2099
years_future_for_model <- if (startsWith(model, "RACMO")) {
  years_future[years_future <= racmo_future_end_year]
} else {
  years_future
}

# The two climatology periods, each tagged with which scenario folder its
# source files live in. `range` always uses the ORIGINAL years_hist/
# years_future vectors for the filename label, even when `years` (what's
# actually used to compute the climatology) has been clipped for a
# particular model.
#
# For the ERA5-driven models we only build a HISTORICAL period, pointed at
# scenario "evaluation" -- the period label and range are kept identical to
# the other models' HISTORICAL entry so output filenames still read
# "..._HISTORICAL_1985_2010.tif" the same way. FUTURE is only computed for
# the non-ERA5 models.
periods <- if (is_era5) {
  list(
    HISTORICAL = list(years = years_hist, scenario = "evaluation",
                      range = paste(min(years_hist), max(years_hist), sep = "_"))
  )
} else {
  list(
    HISTORICAL = list(years = years_hist, scenario = "historical",
                      range = paste(min(years_hist), max(years_hist), sep = "_")),
    FUTURE     = list(years = years_future_for_model, scenario = ssp_scenario,
                      range = paste(min(years_future), max(years_future), sep = "_"))
  )
}

message("Job ", job_index, " -> Model: ", model, " | Variable: temp | Periods: ",
        paste(names(periods), collapse = ", "))

model_dir <- file.path(base_dir, model)

# Same output location as the other scripts, with an extra "comparison"
# subfolder appended.
outpath <- here("Data/Environmental_predictors/PolarRes26", model, "comparison")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

# RACMO and ERA5-driven grids (i.e. everything that isn't HCLIM) carry
# incorrect/missing CRS + extent metadata in their netCDFs. `needs_gridfix`
# flags every non-HCLIM model, and `crsfix` / `exfix` are the corrected
# rotated-pole CRS and extent applied to every raster read for those models
# -- see the fix applied inside load_variable_series() in section 2, right
# after the raster is read in and before anything else touches it.
needs_gridfix <- !startsWith(model, "HCLIM")

exfix <- ext(c(144, 210, -28.1, 25))
crsfix <- "GEOGCRS[\"Rotated_pole\",
    BASEGEOGCRS[\"unknown\",
        DATUM[\"unnamed\",
            ELLIPSOID[\"Sphere\",6371229,0,
                LENGTHUNIT[\"metre\",1,
                    ID[\"EPSG\",9001]]]],
        PRIMEM[\"Greenwich\",0,
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]]],
    DERIVINGCONVERSION[\"Pole rotation (netCDF CF convention)\",
        METHOD[\"Pole rotation (netCDF CF convention)\"],
        PARAMETER[\"Grid north pole latitude (netCDF CF convention)\",5,
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]],
        PARAMETER[\"Grid north pole longitude (netCDF CF convention)\",20,
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]],
        PARAMETER[\"North pole grid longitude (netCDF CF convention)\",0,
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]]],
    CS[ellipsoidal,2],
        AXIS[\"latitude\",north,
            ORDER[1],
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]],
        AXIS[\"longitude\",east,
            ORDER[2],
            ANGLEUNIT[\"degree\",0.0174532925199433,
                ID[\"EPSG\",9122]]]]"

# Template file used to recover a CRS for HCLIM inputs that come out of
# Scripts 1b with an empty/missing CRS (a known HCLIM quirk). Only ever
# used as a source of CRS metadata -- never as data -- and only after the
# extent of the broken input has been checked against this template's
# extent (see section 5a).
# Downloaded as part of script 1a OR you can manually download from:
# https://cordex.dmi.dk/thredds/fileServer/esg_cordex/PolarRes/ANT-12/HCLIMcom-DMI/CESM2/historical/r11i1p1f1/HCLIM43-ALADIN/v1-r1/day/hurs/v20251130/hurs_ANT-12_CESM2_historical_r11i1p1f1_HCLIMcom-DMI_HCLIM43-ALADIN_v1-r1_day_19860101-19901231.nc

hclim_crs_template_path <- here(
  "Data/PolarRes26/HCLIM_CESM2/historical/r11i1p1f1/HCLIM43-ALADIN/v1-r1/day/hurs/v20251130/",
  "hurs_ANT-12_CESM2_historical_r11i1p1f1_HCLIMcom-DMI_HCLIM43-ALADIN_v1-r1_day_19860101-19901231.nc"
)


# ---- 2. Generic engine: find files, load a continuous daily series -------------

# Folder-based match, the variable name must appear as an exact path component (so "tas" never matches
# "tasmax"), and the scenario folder must match.

# This function finds all the files in the folder, matches the ones that are the right variable and scenario (evaluation, historical or ssp) and just returns the matching ones
find_variable_files <- function(model_dir, variable_name, scenario) {
  all_files <- list.files(model_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
  if (length(all_files) == 0) return(character(0))
  components_list <- str_split(all_files, "/")
  is_var  <- vapply(components_list, function(cmp) variable_name %in% cmp, logical(1))
  is_scen <- grepl(scenario, all_files, fixed = TRUE)
  sort(all_files[is_var & is_scen])
}

# Load every matching file into ONE time-ordered SpatRaster with real dates.
# DIAGNOSTICS ASSUMED: no spillover, no gaps, no duplicate dates, daily,
# standard calendar.
load_variable_series <- function(model_dir, variable_name, scenario) {
  files <- find_variable_files(model_dir, variable_name, scenario)
  if (length(files) == 0) {
    stop(sprintf("No files found for variable '%s', scenario '%s' under %s",
                 variable_name, scenario, model_dir))
  }
  r <- terra::rast(files, subds = variable_name)
  
  # RACMO / ERA5-driven grids: fix the CRS + extent metadata in place before
  # anything else (time parsing, unit conversion, etc.) touches this raster.
  if (needs_gridfix) {
    terra::set.crs(r, crsfix)
    terra::set.ext(r, exfix)
  } else if (is.na(terra::crs(r)) || terra::crs(r) == "") {
    
    
    hclim_template <- terra::rast(hclim_crs_template_path)
    
    extents_match <- isTRUE(all.equal(
      as.vector(terra::ext(r)), as.vector(terra::ext(hclim_template)),
      tolerance = 1e-6
    ))
    
    if (!extents_match) {
      stop("Input raster has no CRS and its extent does not match the HCLIM ",
           "CRS template (", hclim_crs_template_path, ") -- refusing to guess ",
           "a CRS for a grid that doesn't line up.",
           "\n  input extent:    ", paste(round(as.vector(terra::ext(r)), 4), collapse = ", "),
           "\n  template extent: ", paste(round(as.vector(terra::ext(hclim_template)), 4), collapse = ", "))
    }
    
    terra::crs(r) <- terra::crs(hclim_template)
    message("  input had no CRS -- extent matched the HCLIM template, so CRS was ",
            "copied from it.")
  }
  
  dates <- as.Date(terra::time(r))
  if (length(dates) != terra::nlyr(r) || anyNA(dates)) {
    stop(sprintf("Could not read a clean daily time dimension for '%s' (%s). ",
                 variable_name, scenario),
         "Re-check Script 1's diagnostics for this variable before trusting this run.")
  }
  ord <- order(dates)
  list(r = r[[ord]], dates = dates[ord])
}

# ---- 3. Generic engine: climatology maths ---------------------------------------

# Per-year mean for one (month, year) combination, averaged across years.
# NOTE: this is the same generic version used for the mean-style variables
# in the full climatology script.
climatological_monthly_mean <- function(r, dates, month, years) {
  yearly <- list()
  for (y in years) {
    idx <- which(lubridate::month(dates) == month & lubridate::year(dates) == y)
    if (length(idx) == 0) {
      warning(sprintf("No data for month %d, year %d -- skipping", month, y))
      next
    }
    sub <- r[[idx]]
    yearly[[as.character(y)]] <- app(sub, mean, na.rm = TRUE)
  }
  if (length(yearly) == 0) stop("No years had any data -- check the date range and scenario coverage.")
  app(rast(yearly), mean, na.rm = TRUE)
}

# All 12 climatological monthly means -- these are each saved individually,
# then used to build the annual and seasonal figures.
climatological_monthly_means_all <- function(r, dates, years) {
  out <- vector("list", 12)
  names(out) <- month.name
  for (m in 1:12) {
    out[[m]] <- climatological_monthly_mean(r, dates, m, years)
  }
  out
}

# Annual figure = mean of the 12 monthly climatologies.
annual_mean_from_monthly <- function(monthly_list) {
  app(rast(monthly_list), mean, na.rm = TRUE)
}

# Mean of a SUBSET of already-computed monthly climatologies -- e.g. DJF or
# JJA. A season's climatology is just the mean of its 3 constituent monthly
# climatologies. Because those monthly climatologies are themselves built
# with no year-crossing (each month is averaged using that same calendar
# year across `years`), this means e.g. "December" here is the
# climatological December for the given years, not the December belonging
# to the following year's Jan/Feb.
seasonal_mean_from_monthly <- function(monthly_list, months) {
  app(rast(monthly_list[months]), mean, na.rm = TRUE)
}

# ---- 4. Unit-conversion helper ---------------------------------------------------
# DIAGNOSTICS ASSUMED: tas is in K, as in the full script.

to_celsius <- function(r) r - 273.15 # K to celsius

# ---- 5. Output writing -----------------------------------------------------------
# All models are written out in their own native grid -- see section 1 for the
# CRS/extent fix applied to every raster as it's read in -- there is no
# resampling onto a common grid.

save_raster <- function(r, filename) {
  path <- file.path(outpath, filename)
  writeRaster(r, path, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  message("  wrote ", filename)
}

# ==============================================================================
# 6. Create the variable
# ==============================================================================

# ---- 1. TEMPERATURE (tas, K -> degC): monthly, annual, summer (DJF), winter (JJA) -
# Run for each period defined in section 1 (HISTORICAL for every model, plus
# FUTURE for the non-ERA5 models).

for (period_name in names(periods)) {
  p <- periods[[period_name]]
  message("-- temp: ", period_name)
  
  series <- load_variable_series(model_dir, "tas", p$scenario) # result of this function is a list of rasters and their dates
  r_celsius <- to_celsius(series$r)
  
  # Annual: 12 monthly climatologies, saved individually, then averaged
  # Inputs are: all the rasters, their dates, and all the years in that period
  monthly <- climatological_monthly_means_all(r_celsius, series$dates, p$years)
  for (m in 1:12) {
    save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Temperature_%s_%s_%s.tif",
                                      month.name[m], period_name, p$range))
  }
  
  # Mean annual = mean of the 12 monthly climatologies
  annual <- annual_mean_from_monthly(monthly)
  save_raster(annual, sprintf("Mean_Annual_Temperature_%s_%s.tif", period_name, p$range))
  
  # Summer (DJF) and winter (JJA) -- both just the mean of the relevant
  # already-computed monthly climatologies, no year-crossing.
  summer <- seasonal_mean_from_monthly(monthly, c("December", "January", "February"))
  save_raster(summer, sprintf("Mean_Summer_Temperature_%s_%s.tif", period_name, p$range))
  
  winter <- seasonal_mean_from_monthly(monthly, c("June", "July", "August"))
  save_raster(winter, sprintf("Mean_Winter_Temperature_%s_%s.tif", period_name, p$range))
}

message("Done: ", model, " / temp / ", paste(names(periods), collapse = ", "))