# ==============================================================================
# PolarRes26 -- STEP 2 (EVALUATION SUBSET): TEMPERATURE ONLY, 1985-2010 & 2071-2100
# ==============================================================================
# Trimmed from the full climatology script: computes ONLY mean-annual
# temperature (plus its 12 monthly climatologies, and DJF/JJA seasonal
# means), across all 7 models, for two periods:
#   - HISTORICAL: 1985-2010, scenario folder depends on the model --
#     "evaluation" for the ERA5-driven models (HCLIM_ERA5, RACMO_ERA5,
#     MetUM_ERA5), "historical" for the rest.
#   - FUTURE: 2071-2100, scenario "ssp370" -- only computed for the
#     non-ERA5 models, since the ERA5-driven models have no ssp series
#     (matching the original full climatology script's logic).
# Output filenames use "HISTORICAL" / "FUTURE" as the period label, and go
# to the usual Data/Environmental_predictors/PolarRes26/<model> location,
# with an extra "comparison" subfolder appended.

# ---- 0. Setup ------------------------------------------------------------------

lib_loc <- paste(getwd(), "/r_lib_new", sep = "")
# .libPaths(lib_loc)

library(terra)
library(here)
library(lubridate)
library(stringr)
library(dplyr, lib.loc = lib_loc)
library(purrr, lib.loc = lib_loc)
library(ncdf4, lib.loc = lib_loc)

# ---- 1. Configuration -----------------------------------------------------------

base_dir <- "/mnt/hpccs01/home/patterc2/n11222026/AntarcticFutureHabitat/Data/PolarRes26"

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5",
            "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5",
            "MetUM_ERA5")

ssp_scenario <- "ssp370"

years_hist   <- seq(1985, 2010, by = 1)
years_future <- seq(2071, 2100, by = 1)

# One job = one model.
args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])

model <- models[job_index]

# ERA5-driven models (HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5) don't have an ssp
# series and their historical-equivalent files live under "evaluation"
# rather than "historical" -- same logic as the original full script.
is_era5 <- endsWith(model, "_ERA5")

# RACMO's ssp series currently ends 2099-12-31 (one year short). For RACMO,
# the years actually used to compute the FUTURE period are clipped to
# 2099, but the output filename keeps the original "2071_2100" label so
# RACMO and HCLIM outputs stay directly comparable.
racmo_future_end_year <- 2099
years_future_for_model <- if (startsWith(model, "RACMO")) {
  years_future[years_future <= racmo_future_end_year]
} else {
  years_future
}

# The two periods this script computes, each tagged with which scenario
# folder its source files live in. FUTURE is only computed for non-ERA5
# models, matching the original full script (ERA5-driven models have no
# ssp series at all).
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
# -- see the fix applied inside load_variable_series(), right after the
# raster is read in and before anything else touches it.
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

# Known HCLIM quirk (see Script 3): some HCLIM outputs come out of Script 1
# with no CRS attached, even though their grid/extent is fine. Only ever
# used as a source of CRS metadata -- never as data -- and only after the
# extent of the broken input has been checked against this template's
# extent (see the fix applied inside load_variable_series()).
hclim_crs_template_path <- here(
  "Data/PolarRes26/HCLIM_CESM2/historical/r11i1p1f1/HCLIM43-ALADIN/v1-r1/day/hurs/v20251130/",
  "hurs_ANT-12_CESM2_historical_r11i1p1f1_HCLIMcom-DMI_HCLIM43-ALADIN_v1-r1_day_19860101-19901231.nc"
)

# ---- 2. Generic engine: find files, load a continuous daily series -------------

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
  
  if (needs_gridfix) {
    terra::set.crs(r, crsfix)
    terra::set.ext(r, exfix)
  } else if (is.na(terra::crs(r)) || terra::crs(r) == "") {
    
    # Known HCLIM quirk: some HCLIM outputs come out of Script 1 with no
    # CRS attached, even though their grid/extent is fine. Recover the CRS
    # from a known-good HCLIM template file -- but only after confirming
    # the extents actually line up, so we never silently stamp a wrong CRS
    # onto a raster that doesn't actually match the template grid.
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

# All 12 climatological monthly means.
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
# JJA. Because the monthly climatologies are built with no year-crossing
# (each month is averaged using that same calendar year across `years`),
# "December" here is the climatological December for the given years, not
# the December belonging to the following year's Jan/Feb.
seasonal_mean_from_monthly <- function(monthly_list, months) {
  app(rast(monthly_list[months]), mean, na.rm = TRUE)
}

# ---- 4. Unit-conversion helper ---------------------------------------------------
# DIAGNOSTICS ASSUMED: tas is in K, as in the full script.

to_celsius <- function(r) r - 273.15 # K to celsius

# ---- 5. Output writing -----------------------------------------------------------

save_raster <- function(r, filename) {
  path <- file.path(outpath, filename)
  writeRaster(r, path, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  message("  wrote ", filename)
}

# ==============================================================================
# 6. TEMPERATURE (tas, K -> degC): mean annual, monthly climatology,
#    summer (DJF), winter (JJA) -- for each period defined above
# ==============================================================================

for (period_name in names(periods)) {
  p <- periods[[period_name]]
  message("-- temp: ", period_name)
  
  series <- load_variable_series(model_dir, "tas", p$scenario)
  r_celsius <- to_celsius(series$r)
  
  # 12 monthly climatologies, saved individually
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