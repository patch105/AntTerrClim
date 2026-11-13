
# ==============================================================================
# HARMONISE THE GRIDS OF THE MONTHLY CHELSA CLIMATOLOGIES 
# ==============================================================================

# Takes every .tif the CHELSA climatology script (2b) wrote to CHELSA/comparison (historical +
# every future model) and produces one masked, regridded output per file:
#   - <name>_ICEFREE.tif   -- reprojected/resampled onto the 1km grid,
#                             masked to the ice-free domain
#
# Core method: (build a target grid, project() with
# bilinear, then mask()) but simplified for this use case:
#   - The target grid is built from Ant_extent_grid_1km.tif instead of
#     coast_domain.tif -- this ONLY supplies extent/CRS/resolution (the
#     1km grid), it is never used to mask/crop anything.
#   - There is no coastline masking step and no separate future ice-free
#     domain -- every file (historical or future) is masked with the
#     same ice_free_domain_1km.tif.
#   - No sea-ice branch -- these are all temperature outputs.
#
# One job = one input file, selected via a single job_index (same pattern
# as Scripts 1/2) -- run with no argument (or an out-of-range one) to
# print the full table of every file that will be processed, and how many
# jobs that is.
# ==============================================================================

# ---- 0. Setup ------------------------------------------------------------------

library(dplyr)
library(purrr)
library(terra)
library(here)

# ---- 1. Configuration -----------------------------------------------------------

# Where Script 1's outputs live (historical + all FUTURE_<model> files).
input_dir <- here("Data/CHELSA/comparison")

# Where this script's regridded outputs go.
output_dir <- here("Data/CHELSA/comparison/regridded")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Target-grid template: extent/CRS/resolution ONLY -- never used to mask.
ant_extent_grid_1km_path <- here("Data/Ant_extent_grid_1km.tif")

# Ice-free mask, applied to every output regardless of historical/future.
ice_free_domain_1km_path <- here("Data/ice_free_domain_1km.tif")

ant_extent_grid_1km <- rast(ant_extent_grid_1km_path)
ice_free_domain_1km <- rast(ice_free_domain_1km_path)

# If TRUE, skip an input file entirely when its output already exists.
# Set to FALSE to always reprocess and overwrite.
SKIP_EXISTING <- TRUE

# ---- 2. Build the job table: every .tif produced ----------------------

job_table <- tibble(
  input_path = list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)
) %>%
  mutate(filename = basename(input_path))

if (nrow(job_table) == 0) {
  stop("No .tif files found under ", input_dir, " -- has Script 1 been run yet?")
}

message(nrow(job_table), " output file(s) found in ", input_dir)

# ---- 3. job_index selects ONE input file  ---------

args <- commandArgs(trailingOnly = TRUE)
job_index <- suppressWarnings(as.integer(args[1]))

if (length(args) == 0 || is.na(job_index) || job_index < 1 || job_index > nrow(job_table)) {
  print(job_table, n = Inf)
  message(nrow(job_table), " job(s) total. Pass a job_index between 1 and ",
          nrow(job_table), " to run a single job.")
  quit(save = "no", status = 0)
}

this_job   <- job_table[job_index, ]
input_path <- this_job$input_path
filename   <- this_job$filename
message("Job ", job_index, "/", nrow(job_table), " -> ", filename)

# ---- 4. Output path ---------------------------------------------------------------

stem <- tools::file_path_sans_ext(filename)
icefree_path <- file.path(output_dir, paste0(stem, "_ICEFREE.tif"))

# ---- 4b. Skip-existing check -------------------------------------------------------

if (SKIP_EXISTING && file.exists(icefree_path)) {
  message("  output already exists -- skipping: ", basename(icefree_path))
  quit(save = "no", status = 0)
}

# ---- 5. Project, resample, mask, save ----------------------------------------------

r <- rast(input_path)

if (is.na(crs(r)) || crs(r) == "") {
  stop("Input raster has no CRS: ", input_path,
       " -- reprojecting without a source CRS would silently produce garbage.")
}

# Step 1: Build a target grid from the 1km extent/CRS/resolution template.
# NOTE: this template supplies extent/CRS/resolution ONLY -- it is never
# used to mask or crop anything (unlike coast_domain in Script 2).
target_grid <- rast(extent = ext(ant_extent_grid_1km), crs = crs(ant_extent_grid_1km),
                    resolution = res(ant_extent_grid_1km))

# Sanity check: extent must be a clean integer multiple of the resolution,
# otherwise project()/rast() will silently produce a grid that doesn't
# tile evenly and downstream resample() will misalign.
stopifnot(
  (xmax(target_grid) - xmin(target_grid)) %% res(target_grid) == 0,
  (ymax(target_grid) - ymin(target_grid)) %% res(target_grid) == 0
)

# Step 2: Reproject into the target grid's CRS/resolution with bilinear interpolation.
r <- project(r, target_grid, method = "bilinear")

# Step 3: Mask to the ice-free domain (same treatment for historical and
# every future model -- no separate future ice-free domain here).
r_icefree <- mask(r, ice_free_domain_1km, maskvalue = NA)
writeRaster(r_icefree, icefree_path, overwrite = TRUE)
message("  wrote ", icefree_path)

message("Done: ", filename)