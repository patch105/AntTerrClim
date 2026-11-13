
# ==============================================================================
# HARMONISE THE GRIDS OF THE MONTHLY ANT AIR ICE CLIMATOLOGIES 
# ==============================================================================

# Takes every .tif the AntAirICE climatology script (2b) wrote to
# AntAirICE/comparison and produces one masked, regridded output per file to match AntTerrClim:
#   - <name>_ICEFREE.tif   -- reprojected/resampled onto the 1km grid,
#                             masked to the ice-free domain
#
# Inputs (HISTORICAL 2003-2014 only -- AntAirICE has no future models):
#   - Climatological_Monthly_Mean_Temperature_<Month>_HISTORICAL_2003_2014.tif (x12)
#   - Mean_Annual_Temperature_HISTORICAL_2003_2014.tif
#   - Mean_Summer_Temperature_HISTORICAL_2003_2014.tif   (DJF)
#   - Mean_Winter_Temperature_HISTORICAL_2003_2014.tif   (JJA)
#   => 15 files expected.
#
# Core method: (build a target grid, project() with bilinear, then mask()):
#   - The target grid is built from Ant_extent_grid_1km.tif -- this ONLY
#     supplies extent/CRS/resolution (the 1km grid), it is never used to
#     mask/crop anything.
#   - Every file is masked with the same ice_free_domain_1km.tif.
#   - No sea-ice branch and no future ice-free domain -- these are all
#     historical temperature outputs.
#
# Runs locally: processes ALL input files in a single run (no job index).
# A failure on one file is reported at the end rather than stopping the rest.
# ==============================================================================

# ---- 0. Setup ------------------------------------------------------------------

library(terra)
library(here)

# ---- 1. Configuration -----------------------------------------------------------

# Where the AntAirICE climatology script's outputs live.
input_dir <- here("Data/AntAirICE/comparison")

# Where this script's regridded outputs go.
output_dir <- here("Data/AntAirICE/comparison/regridded")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Target-grid template: extent/CRS/resolution ONLY -- never used to mask.
ant_extent_grid_1km_path <- here("Data/Ant_extent_grid_1km.tif")

# Ice-free mask, applied to every output.
ice_free_domain_1km_path <- here("Data/ice_free_domain_1km.tif")

ant_extent_grid_1km <- rast(ant_extent_grid_1km_path)
ice_free_domain_1km <- rast(ice_free_domain_1km_path)

# Expected number of AntAirICE climatology outputs (12 monthly + annual +
# summer + winter). Only used for a visibility warning, not a hard stop.
EXPECTED_N_FILES <- 15

# If TRUE, skip an input file when its output already exists.
# Set to FALSE to always reprocess and overwrite.
SKIP_EXISTING <- TRUE

# ---- 2. Find input files ---------------------------------------------------------

# Non-recursive, so the regridded/ subfolder is never picked up as input.
input_paths <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

# Never treat already-regridded outputs as inputs.
input_paths <- input_paths[!grepl("_ICEFREE\\.tif$", input_paths)]

if (length(input_paths) == 0) {
  stop("No .tif files found in ", input_dir,
       " -- has the AntAirICE climatology script been run yet?")
}

message(length(input_paths), " input file(s) found in ", input_dir)

if (length(input_paths) != EXPECTED_N_FILES) {
  warning("Expected ", EXPECTED_N_FILES, " AntAirICE climatology files (12 monthly + ",
          "annual + summer + winter) but found ", length(input_paths), ".")
}

# ---- 3. Build the target grid ONCE -----------------------------------------------

# This template supplies extent/CRS/resolution ONLY -- it is never used to
# mask or crop anything.
target_grid <- rast(extent = ext(ant_extent_grid_1km), crs = crs(ant_extent_grid_1km),
                    resolution = res(ant_extent_grid_1km))

# Sanity check: extent must be a clean integer multiple of the resolution,
# otherwise project()/rast() will silently produce a grid that doesn't
# tile evenly and downstream resample() will misalign.
stopifnot(
  (xmax(target_grid) - xmin(target_grid)) %% res(target_grid) == 0,
  (ymax(target_grid) - ymin(target_grid)) %% res(target_grid) == 0
)

# ---- 4. Process one file ---------------------------------------------------------

process_file <- function(input_path) {
  filename     <- basename(input_path)
  stem         <- tools::file_path_sans_ext(filename)
  icefree_path <- file.path(output_dir, paste0(stem, "_ICEFREE.tif"))
  
  if (SKIP_EXISTING && file.exists(icefree_path)) {
    message("  output already exists -- skipping: ", basename(icefree_path))
    return(invisible("skipped"))
  }
  
  r <- rast(input_path)
  
  if (is.na(crs(r)) || crs(r) == "") {
    stop("Input raster has no CRS: ", input_path,
         " -- reprojecting without a source CRS would silently produce garbage.")
  }
  
  # Step 1: Reproject into the target grid's CRS/resolution with bilinear interpolation.
  r <- project(r, target_grid, method = "bilinear")
  
  # Step 2: Mask to the ice-free domain.
  r_icefree <- mask(r, ice_free_domain_1km, maskvalue = NA)
  writeRaster(r_icefree, icefree_path, overwrite = TRUE)
  message("  wrote ", basename(icefree_path))
  
  invisible("done")
}

# ---- 5. Run all files --------------------------------------------------------------

status <- character(length(input_paths))
names(status) <- basename(input_paths)

for (i in seq_along(input_paths)) {
  message("[", i, "/", length(input_paths), "] ", basename(input_paths[i]))
  status[i] <- tryCatch(
    process_file(input_paths[i]),
    error = function(e) {
      message("  FAILED: ", conditionMessage(e))
      paste("failed:", conditionMessage(e))
    }
  )
  gc()
}

# ---- 6. Summary --------------------------------------------------------------------

message("\nSummary: ",
        sum(status == "done"), " written, ",
        sum(status == "skipped"), " skipped, ",
        sum(startsWith(status, "failed")), " failed.")

if (any(startsWith(status, "failed"))) {
  print(status[startsWith(status, "failed")])
  warning("Some files failed -- see the list above.")
}

message("Done: AntAirICE regrid + ice-free mask")
