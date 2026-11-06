
# ==============================================================================
# SUMMARISE THE 'ADDITIONAL' HCLIM_CESM2 RAW CLIMATE VARIABLES (siconca, rsus, snm) INTO MONTHLY CLIMATOLOGIES
# ==============================================================================

# This script loads the climate variables downloaded with scripts 1a. in a daily format and summarises them to monthly climatologies. 

# The variables are 'additional' ones downloaded separately and of a slightly different structure to most other variables - hence a different summarising method

# Designed to be run on a High-Performance Computer with an array of jobs

# Every job is a separate model and variable - specified from the 'models' and 'variables' lists below

# The script creates a folder directory structured as: here("Data/Environmental_predictors/PolarRes26", model) where the monthly climatologies are saved as .tif files

# The outputs are still in their native extent and resolution. The next step will harmonise them to a common grid.

# * NOTE - each job processes hundreds of rasters and took from 12-40 hours on a high-performance computer


# DIFFERENCES WITH THESE 'ADDITIONAL' FILES:

# Each monthly file's LAST layer is actually the first timestep of the *next*
# month, which has to be dropped before the file's layers can be trusted.

# There are a couple of corrupt files (siconca: Apr/Jun 1994, snc: Nov 1995, rsus: May 1997) which have to be skipped

# rsus (monthly-chunk "additional" files) and rsds (standard continuous
# files) are NOT on exactly the same grid. Before subtracting to get net
# shortwave, rsus is therefore resampled onto the rsds grid (rsds is the
# template) with bilinear interpolation -- see align_to_template() below,
# used by the solar_rad_net and total_solar_rad_net blocks.

# ---- 0. Setup ------------------------------------------------------------------

library(terra)
library(here)
library(lubridate)
library(stringr)

# ---- 1. Configuration -----------------------------------------------------------

model <- "HCLIM_CESM2"  # this script is only ever run for this one model

# Standard PolarRes26 folder (continuous-file variables, e.g. rsds) --
# same base_dir 
base_dir  <- here("Data/PolarRes26")
model_dir <- file.path(base_dir, model)

# The "additional" folder holding siconca / rsus / snm / snc 
# monthly-chunked-file format.
additional_dir <- here("Data/PolarRes26/HCLIM_CESM2_additional")

ssp_scenario <- "ssp370"

years_hist   <- seq(1995, 2014, by = 1)
years_mid    <- seq(2041, 2060, by = 1)
years_future <- seq(2081, 2100, by = 1)

# One job = one variable.
args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])

variables <- list("sea_ice", "solar_rad_net", "mean_melt", "total_melt",
                  "total_solar_rad_net")

variable  <- variables[[job_index]]
message("Job ", job_index, " -> Model: ", model, " | Variable: ", variable)

# Output goes to the SAME place as other HCLIM_CESM2 outputs -- no
# "_additional" in the output path, per the naming convention requested.
outpath <- here("Data/Environmental_predictors/PolarRes26", model)
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

# HCLIM_CESM2 needs no CRS/extent gridfix (that's only for RACMO/ERA5-driven
# grids in the other script), and it isn't ERA5-driven, so it gets the full
# HISTORICAL/MID/FUTURE period set with no RACMO-style future-year clipping.
periods <- list(
  HISTORICAL = list(years = years_hist,   scenario = "historical",
                    range = paste(min(years_hist),   max(years_hist),   sep = "_")),
  MID        = list(years = years_mid,    scenario = ssp_scenario,
                    range = paste(min(years_mid),    max(years_mid),    sep = "_")),
  FUTURE     = list(years = years_future, scenario = ssp_scenario,
                    range = paste(min(years_future), max(years_future), sep = "_"))
)

# ---- 2a. Loader for STANDARD continuous files ----------

find_variable_files <- function(model_dir, variable_name, scenario) {
  all_files <- list.files(model_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
  if (length(all_files) == 0) return(character(0))
  components_list <- str_split(all_files, "/")
  is_var  <- vapply(components_list, function(cmp) variable_name %in% cmp, logical(1))
  is_scen <- grepl(scenario, all_files, fixed = TRUE)
  sort(all_files[is_var & is_scen])
}

load_variable_series <- function(model_dir, variable_name, scenario) {
  files <- find_variable_files(model_dir, variable_name, scenario)
  if (length(files) == 0) {
    stop(sprintf("No files found for variable '%s', scenario '%s' under %s",
                 variable_name, scenario, model_dir))
  }
  r <- rast(files, subds = variable_name)
  dates <- as.Date(terra::time(r))
  if (length(dates) != terra::nlyr(r) || anyNA(dates)) {
    stop(sprintf("Could not read a clean daily time dimension for '%s' (%s).",
                 variable_name, scenario))
  }
  ord <- order(dates)
  list(r = r[[ord]], dates = dates[ord])
}

# ---- 2b. Loader for the "additional" MONTHLY-CHUNK files (siconca/rsus/snm) --

# Matches files by a leading "variable_" prefix on the filename (e.g.
# "rsus_fp_..."), same substring-on-scenario match Script 2 uses.
find_monthly_files <- function(dir, variable_name, scenario) {
  all_files <- list.files(dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
  if (length(all_files) == 0) return(character(0))
  is_var  <- grepl(paste0("^", variable_name, "_"), basename(all_files))
  is_scen <- grepl(scenario, all_files, fixed = TRUE)
  sort(all_files[is_var & is_scen])
}

# See the header comment block at the top of this script for the reasoning
# behind this function. In short: read each monthly file's real per-layer
# dates, keep only the layers that fall in that file's own calendar month
# then concatenate everything into one continuous, date-sorted series.
load_variable_series_monthly_chunks <- function(dir, variable_name, scenario) {
  
  files <- find_monthly_files(dir, variable_name, scenario)
  if (length(files) == 0) {
    stop(sprintf("No monthly-chunk files found for variable '%s', scenario '%s' under %s",
                 variable_name, scenario, dir))
  }
  
  all_r     <- list()
  all_dates <- as.Date(character(0))
  
  for (f in files) {
    
    layers <- tryCatch({
      
      r <- rast(f)
      dates <- as.Date(terra::time(r))
      
      if (length(dates) != terra::nlyr(r) || anyNA(dates)) {
        stop("unreadable/inconsistent time dimension")
      }
      
      # This file's "own" calendar month = the month of its FIRST layer. Any
      # trailing layer(s) dated into the following month are the spillover
      # and get dropped here.
      own_year  <- lubridate::year(dates[1])
      own_month <- lubridate::month(dates[1])
      keep <- lubridate::year(dates) == own_year & lubridate::month(dates) == own_month
      
      if (!any(keep)) stop("no layers within its own month")
      
      list(r = r[[which(keep)]], dates = dates[keep])
      
    }, error = function(e) {
      warning(sprintf("Skipping corrupt/unreadable file (%s):\n  %s", conditionMessage(e), f))
      NULL
    })
    
    if (is.null(layers)) next
    
    all_r[[length(all_r) + 1]] <- layers$r
    all_dates <- c(all_dates, layers$dates)
  }
  
  if (length(all_r) == 0) {
    stop(sprintf("No usable layers assembled for variable '%s', scenario '%s'.",
                 variable_name, scenario))
  }
  
  r_all <- rast(all_r)
  
  # De-duplicate on date (keep first occurrence) in case two files' own-month
  # windows ever overlap, then sort chronologically.
  dup <- duplicated(all_dates)
  if (any(dup)) {
    warning(sprintf("%d duplicate date(s) found for '%s' / '%s' -- keeping first occurrence.",
                    sum(dup), variable_name, scenario))
    r_all     <- r_all[[!dup]]
    all_dates <- all_dates[!dup]
  }
  
  ord <- order(all_dates)
  list(r = r_all[[ord]], dates = all_dates[ord])
}

# Prints which (year, month) combinations are absent from an assembled
# series, for the years actually needed by this period.
# files (siconca: Apr/Jun 1994; snc: Nov 1995; rsus: May 1997) -- anything
# it reports here is a month climatological_monthly_mean()/_sum() will skip
# automatically.

report_missing_year_months <- function(dates, years, variable_name, scenario) {
  present  <- unique(format(dates, "%Y-%m"))
  expected <- as.vector(outer(years, sprintf("%02d", 1:12), paste, sep = "-"))
  missing  <- setdiff(expected, present)
  if (length(missing) > 0) {
    message(sprintf("  [%s / %s] missing %d year-month(s) (auto-skipped by the climatology functions): %s",
                    variable_name, scenario, length(missing), paste(sort(missing), collapse = ", ")))
  } else {
    message(sprintf("  [%s / %s] no missing year-months detected for %d-%d.",
                    variable_name, scenario, min(years), max(years)))
  }
}

# ---- 2c. Grid alignment helper  -----------------------------------------------
# rsus (monthly-chunk files) and rsds (standard continuous files) are not on
# exactly the same grid, so rsus is resampled onto the rsds grid (rsds =
# template) before the two are subtracted to give net shortwave.
#
# This is equivalent to resample(rsus, rsds, method = "bilinear"). Only the
# FIRST layer of rsds is passed as the template because resample() only uses
# the template's geometry (extent, resolution, CRS), which is much cheaper
# than handing it the whole multi-year stack. If the grids already match
# exactly the resample is skipped.

align_to_template <- function(r, template, label = "layer") {
  if (isTRUE(compareGeom(r, template, stopOnError = FALSE))) {
    message("  ", label, " already on the template grid -- no resampling needed.")
    return(r)
  }
  message("  resampling ", label, " onto the rsds grid (bilinear)")
  resample(r, template, method = "bilinear")
}

# ---- 3. Climatology maths  --------------------------------

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

annual_total_from_monthly <- function(monthly_list) {
  app(rast(monthly_list), sum, na.rm = TRUE)
}

seasonal_mean_from_monthly <- function(monthly_list, months) {
  app(rast(monthly_list[months]), mean, na.rm = TRUE)
}

seasonal_total_from_monthly <- function(monthly_list, months) {
  app(rast(monthly_list[months]), sum, na.rm = TRUE)
}

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

climatological_monthly_sums_all <- function(r, dates, years, transform = NULL) {
  out <- vector("list", 12)
  names(out) <- month.name
  for (m in 1:12) {
    out[[m]] <- climatological_monthly_sum(r, dates, m, years, transform = transform)
  }
  out
}

# ---- 4. Unit-conversion helpers (identical to other script) --------------------------

flux_to_mm <- function(r) r * 86400  # daily flux (kg m-2 s-1) to daily mm
to_MJ      <- function(r) r * 0.0864 # W m-2 -> MJ m-2 day-1 (multiply by 86400 s/day, then 1e-6 J->MJ)

# ---- 5. Output writing (identical to other script) ------------------------------------

save_raster <- function(r, filename) {
  path <- file.path(outpath, filename)
  writeRaster(r, path, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  message("  wrote ", filename)
}

# ==============================================================================
# 6. Create the variables
# ==============================================================================

# ---- SEA ICE (siconca, %): monthly climatology, full grid, no masking ----------

if (variable == "sea_ice") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- sea_ice: ", period_name)
    
    series <- load_variable_series_monthly_chunks(additional_dir, "siconca", p$scenario)
    report_missing_year_months(series$dates, p$years, "siconca", p$scenario)
    
    monthly <- climatological_monthly_means_all(series$r, series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Sea_Ice_Concentration_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
  }
}

# ---- SOLAR RADIATION NET, MEAN (net SW = rsds - rsus, W m-2 -> MJ m-2 day-1) -----
# rsds loads the normal (continuous-file) way, exactly as in Script 2.
# rsus loads via the monthly-chunk loader, which will also transparently
# absorb the known corrupt May-1997 rsus file (see report_missing_year_
# months() output below) the same way Script 1's manual NA-fill-then-skip
# did -- that month is simply excluded from the days used to compute the
# climatology for years containing it.
#
# rsus is on a slightly different grid to rsds, so after the two series are
# matched on common dates, rsus is resampled onto the rsds grid (bilinear,
# rsds = template) BEFORE the subtraction.
#
# The MJ conversion (x 0.0864) is applied after differencing rsds - rsus;
# since it's just a scalar multiply, doing it before or after the
# subtraction gives the same result.
if (variable == "solar_rad_net") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- solar_rad_net: ", period_name)
    
    rsds <- load_variable_series(model_dir, "rsds", p$scenario)
    rsus <- load_variable_series_monthly_chunks(additional_dir, "rsus", p$scenario)
    report_missing_year_months(rsus$dates, p$years, "rsus", p$scenario)
    
    common_dates <- as.Date(intersect(rsds$dates, rsus$dates), origin = "1970-01-01")
    if (length(common_dates) == 0) {
      stop("rsds and rsus share no common dates for ", model, " / ", p$scenario)
    }
    
    rsds_r <- rsds$r[[match(common_dates, rsds$dates)]]
    rsus_r <- rsus$r[[match(common_dates, rsus$dates)]]
    
    # Put rsus on the rsds grid before subtracting.
    rsus_r <- align_to_template(rsus_r, rsds_r[[1]], "rsus")
    
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

# ---- SOLAR RADIATION NET, TOTAL (net SW = rsds - rsus, W m-2 -> MJ m-2 day-1) ----
# Same rsds/rsus loading, date-alignment and grid-alignment handling as
# solar_rad_net above (rsds via the standard continuous loader, rsus via the
# monthly-chunk loader with its known May-1997 gap auto-skipped, then rsus
# resampled onto the rsds grid before subtracting). The MJ conversion is
# then applied as the `transform` inside climatological_monthly_sums_all(),
# same method other script uses for total_precip / total_solar_rad_net.
if (variable == "total_solar_rad_net") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_solar_rad_net: ", period_name)
    
    rsds <- load_variable_series(model_dir, "rsds", p$scenario)
    rsus <- load_variable_series_monthly_chunks(additional_dir, "rsus", p$scenario)
    report_missing_year_months(rsus$dates, p$years, "rsus", p$scenario)
    
    common_dates <- as.Date(intersect(rsds$dates, rsus$dates), origin = "1970-01-01")
    if (length(common_dates) == 0) {
      stop("rsds and rsus share no common dates for ", model, " / ", p$scenario)
    }
    
    rsds_r <- rsds$r[[match(common_dates, rsds$dates)]]
    rsus_r <- rsus$r[[match(common_dates, rsus$dates)]]
    
    # Put rsus on the rsds grid before subtracting.
    rsus_r <- align_to_template(rsus_r, rsds_r[[1]], "rsus")
    
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

# ---- MEAN & TOTAL SNOW MELT (snm, converted to mm like precipitation) ------------

if (variable == "mean_melt") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- mean_melt: ", period_name)
    
    series <- load_variable_series_monthly_chunks(additional_dir, "snm", p$scenario)
    report_missing_year_months(series$dates, p$years, "snm", p$scenario)
    r_mm <- flux_to_mm(series$r)
    
    monthly <- climatological_monthly_means_all(r_mm, series$dates, p$years)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Melt_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    annual <- annual_mean_from_monthly(monthly)
    save_raster(annual, sprintf("Mean_Annual_Melt_%s_%s.tif", period_name, p$range))
  }
}

if (variable == "total_melt") {
  
  for (period_name in names(periods)) {
    p <- periods[[period_name]]
    message("-- total_melt: ", period_name)
    
    series <- load_variable_series_monthly_chunks(additional_dir, "snm", p$scenario)
    report_missing_year_months(series$dates, p$years, "snm", p$scenario)
    
    monthly <- climatological_monthly_sums_all(series$r, series$dates, p$years, transform = flux_to_mm)
    for (m in 1:12) {
      save_raster(monthly[[m]], sprintf("Climatological_Monthly_Total_Melt_%s_%s_%s.tif",
                                        month.name[m], period_name, p$range))
    }
    total <- annual_total_from_monthly(monthly)
    save_raster(total, sprintf("Mean_Total_Annual_Melt_%s_%s.tif", period_name, p$range))
  }
}


message("Done: ", model, " / ", variable)