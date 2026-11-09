# ==============================================================================
# SUMMARISE THE RAW CLIMATE VARIABLES INTO MONTHLY CLIMATOLOGIES
# ==============================================================================

# This script loads the climate variables downloaded with scripts 1a. in a daily format and summarises them to monthly climatologies. 

# Designed to be run on a High-Performance Computer with an array of jobs

# Every job is a separate model and variable - specified from the 'models' and 'variables' lists below

# A separate script runs the bioclimatic variables using outputs from this script

# The script creates a folder directory structured as: here("Data/Environmental_predictors/PolarRes26", model) where the monthly climatologies are saved as .tif files

# The outputs are still in their native extent and resolution. The next step will harmonise them to a common grid.

# * NOTE - each job processes hundreds of rasters and took from 12-40 hours on a high-performance computer

# * Also NOTE - a couple of these variables will fail for HCLIM_CESM2,  HCLIM_MPI_ESM1 and HCLIM_ERA5 (ones involving siconca, rsus, snm). That's because they were downloaded separately (in the _Additional scripts) and will be summarised in separate scripts due to slightly different data structures (1b_Summarising_Climate_Vars_HCLIM_CESM2_Additional.R etc.)

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

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5", "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5", "MetUM_ERA5")

ssp_scenario <- "ssp370"

years_hist   <- seq(1995, 2014, by = 1)
years_mid    <- seq(2041, 2060, by = 1)
years_future <- seq(2081, 2100, by = 1)

# One job = one (model, variable) combination, matching the two-argument
# job-array style. args[1] = model index (1-4), args[2] = variable index.
args <- commandArgs(trailingOnly = TRUE)
job_index    <- as.integer(args[1])

variables <- list("temp", "total_DD_minus5", "total_DD_0", "wind", "sea_ice",
                  "total_precip", "total_summer_precip", "mean_precip",
                  "mean_summer_precip", "solar_rad_down", "solar_rad_net","mean_melt", "total_melt", "mean_hurs", "vpd", "temp_min", "temp_max",
                  "total_solar_rad_down", "total_solar_rad_net")

job_table <- expand.grid(model = models, variable = variables,
                         stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

model    <- job_table$model[job_index]
variable <- job_table$variable[job_index]
message("Job ", job_index, " -> Model: ", model, " | Variable: ", variable)

model_dir <- file.path(base_dir, model)

outpath <- here("Data/Environmental_predictors/PolarRes26", model)
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

# RACMO's ssp series currently ends 2099-12-31 (one year short of HCLIM's
# 2100). For RACMO, the years actually used
# to compute the FUTURE period are clipped to 2099, but the output
# filenames keep the original "2081_2100" label (set below) so RACMO and
# HCLIM outputs stay directly comparable.
racmo_future_end_year <- 2099
years_future_for_model <- if (startsWith(model, "RACMO")) {
  years_future[years_future <= racmo_future_end_year]
} else {
  years_future
}

# The three climatology periods, each tagged with which scenario folder its
# source files live in. `range` always uses the ORIGINAL years_hist/mid/
# future vectors for the filename label, even when `years` (what's actually
# used to compute the climatology) has been clipped for a particular model.
#
# ERA5-driven models (HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5) don't have an ssp
# series at all, and their historical-equivalent files live under an
# "evaluation" folder rather than "historical". So for these models we only
# build a HISTORICAL period, pointed at scenario "evaluation" -- the period
# label and range are kept identical to the other models' HISTORICAL entry
# so output filenames still read "..._HISTORICAL_1995_2014.tif" the same way.
is_era5 <- endsWith(model, "_ERA5")

periods <- if (is_era5) {
  list(
    HISTORICAL = list(years = years_hist, scenario = "evaluation",
                      range = paste(min(years_hist), max(years_hist), sep = "_"))
  )
} else {
  list(
    HISTORICAL = list(years = years_hist,   scenario = "historical",
                      range = paste(min(years_hist), max(years_hist), sep = "_")),
    MID        = list(years = years_mid,    scenario = ssp_scenario,
                      range = paste(min(years_mid), max(years_mid), sep = "_")),
    FUTURE     = list(years = years_future_for_model, scenario = ssp_scenario,
                      range = paste(min(years_future), max(years_future), sep = "_"))
  )
}

# ---- 2. Generic engine: find files, load a continuous daily series -------------

# Folder-based match, the variable name must appear as an exact path component (so "tas" never matches
# "tasmax"), and the scenario folder must match.

# This function finds all the files in the folder, matches the ones that are the right variable and scenario (historical or ssp) and just returns the matching ones
find_variable_files <- function(model_dir, variable_name, scenario) {
  
  all_files <- list.files(model_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
  
  if (length(all_files) == 0) return(character(0))
  components_list <- str_split(all_files, "/")
  
  is_var   <- vapply(components_list, function(cmp) variable_name %in% cmp, logical(1))
  
  is_scen  <- grepl(scenario, all_files, fixed = TRUE)
  
  sort(all_files[is_var & is_scen])
}

# Load every matching file into ONE time-ordered SpatRaster with real dates.

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
# NOTE: this is the generic version used by all monthly-climatology
# variables, including sea ice, temp, wind, precip, solar rad, melt, snow,
# hurs and vpd.
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

# All 12 climatological monthly means -- used for the "mean annual X"
# variables, which also save each month.
climatological_monthly_means_all <- function(r, dates, years) {
  out <- vector("list", 12)
  names(out) <- month.name
  for (m in 1:12) {
    out[[m]] <- climatological_monthly_mean(r, dates, m, years)
  }
  out
}

annual_mean_from_monthly <- function(monthly_list) {
  app(rast(monthly_list), mean, na.rm = TRUE)
}

# Sum (not mean) of a 12-month climatology list -- used for the "total"
# variables (total_DD, total_precip, total_melt, total_solar_rad_down,
# total_solar_rad_net) so their annual figure is a true Jan-Dec sum of the
# monthly climatological totals, rather than an average of them.
annual_total_from_monthly <- function(monthly_list) {
  app(rast(monthly_list), sum, na.rm = TRUE)
}

# Mean of a SUBSET of already-computed monthly climatologies -- e.g. DJF or
# JJA. This is the general building block for every seasonal "mean" style
# output (temp, mean_hurs, mean_summer_precip): a season's climatology is
# just the mean of its 3 constituent monthly climatologies. Because those
# monthly climatologies are themselves built with no year-crossing (each
# month is averaged using that same calendar year across `years`), this
# means e.g. "December" here is the climatological December for the given
# years, not the December belonging to the following year's Jan/Feb.
seasonal_mean_from_monthly <- function(monthly_list, months) {
  app(rast(monthly_list[months]), mean, na.rm = TRUE)
}

# Sum of a SUBSET of already-computed monthly climatology TOTALS -- the
# "total" counterpart to seasonal_mean_from_monthly(), used for
# total_summer_precip.
seasonal_total_from_monthly <- function(monthly_list, months) {
  app(rast(monthly_list[months]), sum, na.rm = TRUE)
}

# Per-year monthly TOTAL for one (month, year) combination, averaged
# across years -- the "total" counterpart to climatological_monthly_mean()
# above. This is what lets the accumulation/total-style variables
# (degree-days, melt, precip, solar rad) be built the same way as the
# mean-style ones: 12 climatological monthly figures (here, monthly SUMS
# rather than monthly MEANS), saved individually, then combined with
# annual_total_from_monthly() for the annual figure.
climatological_monthly_sum <- function(r, dates, month, years, transform = NULL) {
  yearly <- list()
  for (y in years) {
    idx <- which(lubridate::month(dates) == month & lubridate::year(dates) == y)
    if (length(idx) == 0) {
      warning(sprintf("No data for month %d, year %d -- skipping", month, y))
      next
    }
    sub <- r[[idx]]
    if (!is.null(transform)) sub <- transform(sub)
    yearly[[as.character(y)]] <- app(sub, sum, na.rm = TRUE)
  }
  if (length(yearly) == 0) stop("No years had any data -- check the date range and scenario coverage.")
  app(rast(yearly), mean, na.rm = TRUE)
}

# All 12 climatological monthly TOTALS -- mirrors
# climatological_monthly_means_all() but summing within each month instead
# of averaging.
climatological_monthly_sums_all <- function(r, dates, years, transform = NULL) {
  out <- vector("list", 12)
  names(out) <- month.name
  for (m in 1:12) {
    out[[m]] <- climatological_monthly_sum(r, dates, m, years, transform = transform)
  }
  out
}

# ---- 4. Unit-conversion helpers --------------------------------------------------

to_celsius <- function(r) r - 273.15 # K to celcius
flux_to_mm <- function(r) r * 86400 # Daily flux to daily mm by multiplying by seconds
to_MJ <- function(r) r * 0.0864 # W m-2 -> MJ m-2 day-1 (multiply by 86400 s/day, then 1e-6 J->MJ)

# ---- 5. Output writing -----------------------------------------------------------
# RACMO and ERA5-driven models are written out in their own (fixed) native grid -- see
# section 1 for the CRS/extent fix applied to every raster as it's read in
# -- there is no resampling onto a common HCLIM grid.

save_raster <- function(r, filename) {
  path <- file.path(outpath, filename)
  writeRaster(r, path, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  message("  wrote ", filename)
}

# ==============================================================================
# 7. Create the variables 
# ==============================================================================

# ---- 1. TEMPERATURE (tas, K -> degC): annual, summer (DJF), winter (JJA) -------
if (variable == "temp") {
  
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
    annual <- annual_mean_from_monthly(monthly)
    save_raster(annual, sprintf("Mean_Annual_Temperature_%s_%s.tif", period_name, p$range))
    
    # Summer (DJF) and winter (JJA) -- both just the mean of the relevant
    # already-computed monthly climatologies, no year-crossing.
    summer <- seasonal_mean_from_monthly(monthly, c("December", "January", "February"))
    save_raster(summer, sprintf("Mean_Summer_Temperature_%s_%s.tif", period_name, p$range))
    
    winter <- seasonal_mean_from_monthly(monthly, c("June", "July", "August"))
    save_raster(winter, sprintf("Mean_Winter_Temperature_%s_%s.tif", period_name, p$range))
  }
}

# ---- 1b. TASMIN  (K -> degC): monthly climatology only ------------------
# Added for the separate bioclim-summarising script, which needs monthly
# min/max temperature climatologies (not the "mean" tas used by the "temp"
# variable above). 
if (variable == "temp_min") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- temp_min: ", period_name)
    
    tasmin_series <- load_variable_series(model_dir, "tasmin", p$scenario)
    
    tasmin_celsius <- to_celsius(tasmin_series$r)
    
    monthly_min <- climatological_monthly_means_all(tasmin_celsius, tasmin_series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly_min[[m]], sprintf("Climatological_Monthly_Mean_TasMin_%s_%s_%s.tif",
                                            month.name[m], period_name, p$range))
    }
    
  }
}

# ---- 1c. TASMAX (K -> degC): monthly climatology only ------------------
# Added for the separate bioclim-summarising script, which needs monthly
# min/max temperature climatologies (not the "mean" tas used by the "temp"
# variable above).
if (variable == "temp_max") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- temp_max: ", period_name)
    
    tasmax_series <- load_variable_series(model_dir, "tasmax", p$scenario)
    
    tasmax_celsius <- to_celsius(tasmax_series$r)
    
    monthly_min <- climatological_monthly_means_all(tasmax_celsius, tasmax_series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly_min[[m]], sprintf("Climatological_Monthly_Mean_TasMax_%s_%s_%s.tif",
                                            month.name[m], period_name, p$range))
    }
    
  }
}

# ---- 2. TOTAL POSITIVE DEGREE-DAYS, two thresholds: -5 degC and 0 degC (tas) ----
if (variable == "total_DD_minus5") {
  
  dd_transform <- function(r_daily_celsius, limit) {
    ifel(r_daily_celsius > limit, r_daily_celsius - limit, 0)
  }
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_DD_minus5: ", period_name)
    
    series <- load_variable_series(model_dir, "tas", p$scenario)
    r_celsius <- to_celsius(series$r)
    
    monthly_minus5 <- climatological_monthly_sums_all(r_celsius, series$dates, p$years,
                                                      transform = function(r) dd_transform(r, -5))
    for (m in 1:12) {
      save_raster(monthly_minus5[[m]], sprintf("Climatological_Monthly_Total_Degree_Days-5_%s_%s_%s.tif",
                                               month.name[m], period_name, p$range))
    }
    dd_minus5 <- annual_total_from_monthly(monthly_minus5)
    save_raster(dd_minus5, sprintf("Mean_Total_Annual_Degree_Days-5_%s_%s.tif", period_name, p$range))
  }
}

if (variable == "total_DD_0") {
  
  dd_transform <- function(r_daily_celsius, limit) {
    ifel(r_daily_celsius > limit, r_daily_celsius - limit, 0)
  }
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_DD_0: ", period_name)
    
    series <- load_variable_series(model_dir, "tas", p$scenario)
    r_celsius <- to_celsius(series$r)
    
    monthly_0 <- climatological_monthly_sums_all(r_celsius, series$dates, p$years,
                                                 transform = function(r) dd_transform(r, 0))
    for (m in 1:12) {
      save_raster(monthly_0[[m]], sprintf("Climatological_Monthly_Total_Degree_Days0_%s_%s_%s.tif",
                                          month.name[m], period_name, p$range))
    }
    dd_0 <- annual_total_from_monthly(monthly_0)
    save_raster(dd_0, sprintf("Mean_Total_Annual_Degree_Days0_%s_%s.tif", period_name, p$range))
  }
}

# ---- 3. WIND SPEED (sfcWind, m s-1, no conversion): annual ----------------------
if (variable == "wind") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- wind: ", period_name)
    
    series <- load_variable_series(model_dir, "sfcWind", p$scenario)
    
    monthly <- climatological_monthly_means_all(series$r, series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Wind_Speed_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    save_raster(annual, sprintf("Mean_Annual_Wind_Speed_%s_%s.tif", period_name, p$range))
  }
}

# ---- 4. SEA ICE (siconca, %): monthly climatology, full grid -------------------
if (variable == "sea_ice") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- sea_ice: ", period_name)
    
    series <- load_variable_series(model_dir, "siconca", p$scenario)
    
    monthly <- climatological_monthly_means_all(series$r, series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Sea_Ice_Concentration_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
  }
}

# ---- 5a. TOTAL ANNUAL PRECIPITATION (pr, kg m-2 s-1 -> mm) -----------------------
if (variable == "total_precip") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_precip: ", period_name)
    
    series <- load_variable_series(model_dir, "pr", p$scenario)
    
    # daily -> monthly sum -> monthly climatology, saved individually
    monthly <- climatological_monthly_sums_all(series$r, series$dates, p$years, transform = flux_to_mm)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Total_Precipitation_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    total <- annual_total_from_monthly(monthly)
    save_raster(total, sprintf("Total_Annual_Precipitation_%s_%s.tif", period_name, p$range))
  }
}

# ---- 5b. TOTAL SUMMER (DJF) PRECIPITATION (pr, mm) -------------------------------
if (variable == "total_summer_precip") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_summer_precip: ", period_name)
    
    series <- load_variable_series(model_dir, "pr", p$scenario)
    
    monthly_summer <- list(
      December = climatological_monthly_sum(series$r, series$dates, 12, p$years, transform = flux_to_mm),
      January  = climatological_monthly_sum(series$r, series$dates, 1,  p$years, transform = flux_to_mm),
      February = climatological_monthly_sum(series$r, series$dates, 2,  p$years, transform = flux_to_mm)
    )
    total <- seasonal_total_from_monthly(monthly_summer, c("December", "January", "February"))
    save_raster(total, sprintf("Mean_Total_Summer_Precipitation_%s_%s.tif", period_name, p$range))
  }
}

# ---- 5c. MEAN ANNUAL PRECIPITATION (pr, mm/day): annual --------------------------
if (variable == "mean_precip") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- mean_precip: ", period_name)
    
    series <- load_variable_series(model_dir, "pr", p$scenario)
    r_mm <- flux_to_mm(series$r)
    
    monthly <- climatological_monthly_means_all(r_mm, series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Precipitation_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    save_raster(annual, sprintf("Mean_Annual_Precipitation_%s_%s.tif", period_name, p$range))
  }}

# ---- 5d. MEAN SUMMER (DJF) PRECIPITATION (pr, mm/day) ----------------------------
if (variable == "mean_summer_precip") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- mean_summer_precip: ", period_name)
    
    series <- load_variable_series(model_dir, "pr", p$scenario)
    r_mm <- flux_to_mm(series$r)
    
    monthly_summer <- list(
      December = climatological_monthly_mean(r_mm, series$dates, 12, p$years),
      January  = climatological_monthly_mean(r_mm, series$dates, 1,  p$years),
      February = climatological_monthly_mean(r_mm, series$dates, 2,  p$years)
    )
    summer <- seasonal_mean_from_monthly(monthly_summer, c("December", "January", "February"))
    save_raster(summer, sprintf("Mean_Summer_Precipitation_%s_%s.tif", period_name, p$range))
  }
}

# ---- 6a. SOLAR RADIATION NET, MEAN (net SW = rsds - rsus, W m-2 -> MJ m-2 day-1) --
# Converted to MJ m-2 day-1 (x 0.0864) after differencing rsds - rsus; since
# the conversion is just a scalar multiply, doing it before or after the
# subtraction gives the same result.
if (variable == "solar_rad_net") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- solar_rad_net: ", period_name)
    
    rsds <- load_variable_series(model_dir, "rsds", p$scenario)
    rsus <- load_variable_series(model_dir, "rsus", p$scenario)
    
    # rsus has a ~15-year gap (2060-2075) that rsds doesn't have. We don't
    # care about that gap for the climatology (missing years are simply
    # skipped by climatological_monthly_mean), but we do need the two
    # rasters to be subtracted date-for-date. So instead of requiring the
    # full series to line up, restrict both to the dates they have in
    # common before differencing.
    common_dates <- as.Date(intersect(rsds$dates, rsus$dates), origin = "1970-01-01")
    
    if (length(common_dates) == 0) {
      stop("rsds and rsus share no common dates for ", model, " / ", p$scenario)
    }
    
    rsds_r <- rsds$r[[match(common_dates, rsds$dates)]]
    rsus_r <- rsus$r[[match(common_dates, rsus$dates)]]
    
    swnet_r <- to_MJ(rsds_r - rsus_r)
    
    monthly <- climatological_monthly_means_all(swnet_r, common_dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Net_Solar_Radiation_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    save_raster(annual, sprintf("Mean_Annual_Net_Solar_Radiation_%s_%s.tif", period_name, p$range))
  }
}


# ---- 6b. SOLAR RADIATION DOWNWELLING, MEAN (rsds, W m-2 -> MJ m-2 day-1) --------

if (variable == "solar_rad_down") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- solar_rad_down: ", period_name)
    
    rsds <- load_variable_series(model_dir, "rsds", p$scenario)
    rsds_MJ <- to_MJ(rsds$r)
    
    monthly <- climatological_monthly_means_all(rsds_MJ, rsds$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Downwelling_Solar_Radiation_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    save_raster(annual, sprintf("Mean_Annual_Downwelling_Solar_Radiation_%s_%s.tif", period_name, p$range))
  }
}  

# ---- 6c. SOLAR RADIATION DOWNWELLING, TOTAL (rsds, W m-2 -> MJ m-2 day-1) -------

if (variable == "total_solar_rad_down") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_solar_rad_down: ", period_name)
    
    rsds <- load_variable_series(model_dir, "rsds", p$scenario)
    
    monthly <- climatological_monthly_sums_all(rsds$r, rsds$dates, p$years, transform = to_MJ)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Total_Downwelling_Solar_Radiation_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    total <- annual_total_from_monthly(monthly)
    save_raster(total, sprintf("Total_Annual_Downwelling_Solar_Radiation_%s_%s.tif", period_name, p$range))
  }
}

# ---- 6d. SOLAR RADIATION NET, TOTAL (net SW = rsds - rsus, W m-2 -> MJ m-2 day-1) --
# Same date-alignment handling as solar_rad_net (6a) -- rsus has a gap that
# rsds doesn't, so both series are restricted to their common dates before
# differencing. The MJ conversion is then applied as the `transform` inside
# climatological_monthly_sums_all(), same method as total_precip.

if (variable == "total_solar_rad_net") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_solar_rad_net: ", period_name)
    
    rsds <- load_variable_series(model_dir, "rsds", p$scenario)
    rsus <- load_variable_series(model_dir, "rsus", p$scenario)
    
    common_dates <- as.Date(intersect(rsds$dates, rsus$dates), origin = "1970-01-01")
    
    if (length(common_dates) == 0) {
      stop("rsds and rsus share no common dates for ", model, " / ", p$scenario)
    }
    
    rsds_r <- rsds$r[[match(common_dates, rsds$dates)]]
    rsus_r <- rsus$r[[match(common_dates, rsus$dates)]]
    
    swnet_r <- rsds_r - rsus_r
    
    monthly <- climatological_monthly_sums_all(swnet_r, common_dates, p$years, transform = to_MJ)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Total_Net_Solar_Radiation_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    total <- annual_total_from_monthly(monthly)
    save_raster(total, sprintf("Total_Annual_Net_Solar_Radiation_%s_%s.tif", period_name, p$range))
  }
}


# ---- 7. MEAN & TOTAL SNOW MELT (snm, converted to mm like precipitation) -------
# Mass of ice/snow melted at the surface per second (magnitude of melt)
if (variable == "mean_melt") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- mean_melt: ", period_name)
    
    series <- load_variable_series(model_dir, "snm", p$scenario)
    r_mm <- flux_to_mm(series$r)  # now converted to mm, same as precipitation
    
    monthly <- climatological_monthly_means_all(r_mm, series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Melt_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    save_raster(annual, sprintf("Mean_Annual_Melt_%s_%s.tif", period_name, p$range))
  }
}

# ---- 7b. TOTAL ANNUAL SNOW MELT (snm, converted to mm like precipitation) -------
if (variable == "total_melt") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_melt: ", period_name)
    
    series <- load_variable_series(model_dir, "snm", p$scenario)
    
    monthly <- climatological_monthly_sums_all(series$r, series$dates, p$years, transform = flux_to_mm)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Total_Melt_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    total <- annual_total_from_monthly(monthly)
    save_raster(total, sprintf("Mean_Total_Annual_Melt_%s_%s.tif", period_name, p$range))
  }
}


# ---- 8. MEAN RELATIVE HUMIDITY (hurs, %): monthly, annual, summer, winter ------

if (variable == "mean_hurs") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- mean_hurs: ", period_name)
    
    series <- load_variable_series(model_dir, "hurs", p$scenario)
    
    monthly <- climatological_monthly_means_all(series$r, series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Relative_Humidity_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    save_raster(annual, sprintf("Mean_Annual_Relative_Humidity_%s_%s.tif", period_name, p$range))
    
    summer <- seasonal_mean_from_monthly(monthly, c("December", "January", "February"))
    save_raster(summer, sprintf("Mean_Summer_Relative_Humidity_%s_%s.tif", period_name, p$range))
    
    winter <- seasonal_mean_from_monthly(monthly, c("June", "July", "August"))
    save_raster(winter, sprintf("Mean_Winter_Relative_Humidity_%s_%s.tif", period_name, p$range))
  }
}

# ---- 9. VAPOUR PRESSURE DEFICIT (Pa): monthly, annual --------------------------
# Computed from hurs (relative humidity, %) and tas (air temperature), daily,
# then run through the same monthly-climatology + annual-mean pattern as
# temp/wind/etc.
#
# This uses the same Sonntag 1990 / WMO(2008) Magnus formula as
# bigleaf::rH.to.VPD() / Esat.slope() -- but reimplemented directly as raster
# arithmetic below (compute_vpd_Pa()) rather than calling those functions on
# the rasters directly. bigleaf's Esat.slope() ends with
# `data.frame(Esat, Delta)`, which doesn't work on a SpatRaster, so calling
# rH.to.VPD() on raster inputs would error. 
if (variable == "vpd") {
  
  # rH.to.VPD(rH, Tair) returns VPD in kPa; hurs is %, so rH = hurs/100
  compute_vpd_Pa <- function(hurs_pct, tas_celsius) {
    a <- 611.2; b <- 17.62; c <- 243.12  # Sonntag 1990 Magnus coefficients
    rH <- hurs_pct / 100
    esat_Pa <- a * exp((b * tas_celsius) / (c + tas_celsius))
    vpd_Pa <- esat_Pa * (1 - rH)
    vpd_Pa
  }
  
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- vpd: ", period_name)
    
    hurs_series <- load_variable_series(model_dir, "hurs", p$scenario)
    tas_series  <- load_variable_series(model_dir, "tas",  p$scenario)
    
    if (!identical(hurs_series$dates, tas_series$dates)) {
      stop("hurs and tas dates don't line up for ", model, " / ", p$scenario,
           " -- check Script 1's date-gap report for these two variables.")
    }
    
    tas_celsius <- to_celsius(tas_series$r)
    vpd_daily_Pa <- compute_vpd_Pa(hurs_series$r, tas_celsius)
    
    monthly <- climatological_monthly_means_all(vpd_daily_Pa, hurs_series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_VPD_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    save_raster(annual, sprintf("Mean_Annual_VPD_%s_%s.tif", period_name, p$range))
  }
}

message("Done: ", model, " / ", variable)