
# ==============================================================================
# CREAT MONTHLY CLIMATOLOGIES FOR THE ANTAIR ICE DATASET TO COMPARE WITH ANTTERRCLIM
# ==============================================================================

# Takes every daily AntAir ICE .tif file which has been manually downloaded previously and summarises them into monthly climatologies following the methods of script 1b

# It also calculates the mean annual temperature, mean austral summer temperature and mean austral winter temperature

# This was run on a local computer and took about 5 hours

# The outputs are saved in the "Data/AntTerrClim" folder

# Single job and the rasters are loaded per (month, year) so that only the ~28-31 daily files needed for one month/year combination are ever in memory at once.

# No reprojection and no ice-free masking are applied in this script.

# NA HANDLING:
#  - Some AntAirICE daily files are known to potentially contain NA cells. This script:
#     1. Skips unreadable/corrupt files with a warning rather than crashing.
#     2. Computes per-pixel means with na.rm = TRUE, so a day with NA at a
#        given pixel doesn't poison that pixel's month/year mean.
#     3. Explicitly converts any resulting NaN (pixels that were NA on
#        EVERY day available) back to proper NA, since mean(..., na.rm=TRUE)
#        over an all-NA vector returns NaN, not NA, and the two are not
#        always interchangeable downstream (e.g. terra I/O, some GDAL NA
#        handling).
#     4. Reports per (month, year) how many daily files were found/read
#        successfully vs how many days that month should have, and warns
#        if coverage is low, so sparse months are visible rather than silent.

library(terra)
library(here)
library(arrow)
library(lubridate)

# ---- 1. Configuration ------------------------------------------------------

years_eval   <- seq(2003, 2014, by = 1)
period_name  <- "HISTORICAL"
period_range <- paste(min(years_eval), max(years_eval), sep = "_")

# Minimum fraction of a month's calendar days that must have a successfully
# read file before we proceed without extra warning (still computes the mean
# either way -- this is a visibility/QA threshold, not a hard stop).
min_coverage_frac <- 0.5

data_dir <- here("Data/AntAirICE")
outpath  <- here("Data/AntAirICE/comparison")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

# ---- 2. File discovery & date parsing --------------------------------------
# Filenames look like: AntAir_ICE_2003_001.tif  (year, day-of-year)
# Build a lookup table (date, path) -- no rasters read yet, just filenames.

variable_paths <- list.files(
  data_dir,
  pattern    = "\\.tif$",
  full.names = TRUE,
  recursive  = TRUE
)

is_icefree <- grepl("_ICEFREE\\.tif$", variable_paths)
if (any(is_icefree)) {
  message("Excluding ", sum(is_icefree), " already-processed '_ICEFREE.tif' file(s) from the summary.")
}
variable_paths <- variable_paths[!is_icefree]

variable_names <- basename(variable_paths)

year_str <- sub(".*_(\\d{4})_(\\d{3})\\.tif$", "\\1", variable_names)
doy_str  <- sub(".*_(\\d{4})_(\\d{3})\\.tif$", "\\2", variable_names)
years_all <- as.integer(year_str)
doy_all   <- as.integer(doy_str)

if (any(is.na(years_all)) || any(is.na(doy_all))) {
  stop("Could not parse year/day-of-year from one or more AntAirICE filenames -- ",
       "check the naming convention (expected AntAir_ICE_<year>_<doy>.tif).")
}

dates_all <- as.Date(doy_all - 1, origin = paste0(years_all, "-01-01"))

# Restrict to the 2003-2014 evaluation period
keep <- years_all %in% years_eval
file_index <- data.frame(
  path  = variable_paths[keep],
  date  = dates_all[keep],
  year  = years_all[keep],
  month = lubridate::month(dates_all[keep]),
  stringsAsFactors = FALSE
)
file_index <- file_index[order(file_index$date), ]

message("Found ", nrow(file_index), " AntAirICE daily rasters for ", period_range, ".")

# ---- 3. Helpers for safe loading & NA handling ------------------------------

# Safely read layer 1 of a single file. Returns NULL (with a warning) if the
# file can't be read at all, rather than crashing the whole run.
safe_read_layer1 <- function(path) {
  out <- tryCatch(
    terra::rast(path)[[1]],
    error = function(e) {
      warning(sprintf("Could not read '%s' -- skipping this file (%s)",
                      path, conditionMessage(e)))
      NULL
    }
  )
  out
}

# Convert NaN -> NA in a SpatRaster. app()'s mean(..., na.rm=TRUE) returns
# NaN for pixels that were NA in every input layer; this puts a proper NA
# back so downstream tools treat "no data" consistently.
nan_to_na <- function(r) {
  terra::ifel(is.nan(r), NA, r)
}

# ---- 4. Climatology engine (loads per month/year only) ---------------------

# Per-year mean for one (month, year) combination -- loads ONLY the files
# belonging to that month/year, computes the mean, then those rasters can
# be garbage-collected before moving to the next year.
climatological_monthly_mean <- function(file_index, month, years, min_coverage_frac) {
  yearly <- list()
  
  for (y in years) {
    sub_idx <- file_index[file_index$month == month & file_index$year == y, ]
    
    if (nrow(sub_idx) == 0) {
      warning(sprintf("No files found for month %d, year %d -- skipping this year", month, y))
      next
    }
    
    # Safely load only this month/year's daily files, layer 1 each,
    # dropping any that fail to read.
    layers <- lapply(sub_idx$path, safe_read_layer1)
    layers <- layers[!vapply(layers, is.null, logical(1))]
    
    if (length(layers) == 0) {
      warning(sprintf("All files for month %d, year %d failed to read -- skipping this year", month, y))
      next
    }
    
    # Coverage check: how many days of that calendar month were we able to
    # read, vs how many the month should have (for this year).
    days_in_this_month <- lubridate::days_in_month(as.Date(sprintf("%04d-%02d-01", y, month)))
    coverage_frac <- length(layers) / days_in_this_month
    if (coverage_frac < min_coverage_frac) {
      warning(sprintf(
        "Low coverage for month %d, year %d: only %d/%d days readable (%.0f%%)",
        month, y, length(layers), days_in_this_month, coverage_frac * 100
      ))
    }
    
    r_month <- rast(layers)
    
    # Scale: multiply by 0.1 (same unit conversion as the daily script)
    r_month <- r_month * 0.1
    
    # Per-pixel mean across the days available this month/year, ignoring
    # per-pixel NAs. Pixels NA on every available day become NaN here --
    # fixed immediately below.
    r_month_mean <- app(r_month, mean, na.rm = TRUE)
    r_month_mean <- nan_to_na(r_month_mean)
    
    yearly[[as.character(y)]] <- r_month_mean
    
    rm(r_month, layers); gc()
  }
  
  if (length(yearly) == 0) stop("No years had any usable data -- check date range/file coverage.")
  
  out <- app(rast(yearly), mean, na.rm = TRUE)
  nan_to_na(out)
}

# All 12 climatological monthly means.
climatological_monthly_means_all <- function(file_index, years, min_coverage_frac) {
  out <- vector("list", 12)
  names(out) <- month.name
  for (m in 1:12) {
    message("  computing climatology for ", month.name[m], " ...")
    out[[m]] <- climatological_monthly_mean(file_index, m, years, min_coverage_frac)
  }
  out
}

# Annual figure = mean of the 12 monthly climatologies.
annual_mean_from_monthly <- function(monthly_list) {
  out <- app(rast(monthly_list), mean, na.rm = TRUE)
  nan_to_na(out)
}

# Mean of a subset of already-computed monthly climatologies (DJF or JJA).
# No year-crossing: "December" is the climatological December for
# years_eval, not the Dec belonging to the following year's Jan/Feb.
seasonal_mean_from_monthly <- function(monthly_list, months) {
  out <- app(rast(monthly_list[months]), mean, na.rm = TRUE)
  nan_to_na(out)
}

# ---- 5. Output writing -------------------------------------------------------

save_raster <- function(r, filename) {
  path <- file.path(outpath, filename)
  writeRaster(r, path, gdal = c("COMPRESS=NONE"), overwrite = TRUE)
  message("  wrote ", filename)
}

# ---- 6. Run: monthly climatologies, mean annual, summer (DJF), winter (JJA) --

message("-- AntAirICE temp: ", period_name)

monthly <- climatological_monthly_means_all(file_index, years_eval, min_coverage_frac)
for (m in 1:12) {
  save_raster(monthly[[m]], sprintf("Climatological_Monthly_Mean_Temperature_%s_%s_%s.tif",
                                    month.name[m], period_name, period_range))
}

annual <- annual_mean_from_monthly(monthly)
save_raster(annual, sprintf("Mean_Annual_Temperature_%s_%s.tif", period_name, period_range))

summer <- seasonal_mean_from_monthly(monthly, c("December", "January", "February"))
save_raster(summer, sprintf("Mean_Summer_Temperature_%s_%s.tif", period_name, period_range))

winter <- seasonal_mean_from_monthly(monthly, c("June", "July", "August"))
save_raster(winter, sprintf("Mean_Winter_Temperature_%s_%s.tif", period_name, period_range))

message("Done: AntAirICE / temp / ", period_name)

