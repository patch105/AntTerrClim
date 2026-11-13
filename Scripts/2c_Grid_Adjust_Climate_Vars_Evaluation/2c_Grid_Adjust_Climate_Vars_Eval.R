
# ==============================================================================
# HARMONISE THE GRIDS OF THE MONTHLY CLIMATOLOGIES 
# ==============================================================================

# Same grid-adjustment logic as the full 1b script, but scoped to the
# "comparison" outputs (temperature mean annual, monthly
# climatologies, summer, winter) produced by the 2b script. Standardises the climatologies to a common grid of 1 km resolution for comparison with Ant Air ICE and CHELSA.
# Reads from each model's Data/PolarRes26/<model>/comparison folder, and
# writes to a matching .../Regridded/<model>/comparison folder.

# ---- 0. Setup ------------------------------------------------------------------

library(dplyr)
library(purrr)
library(terra)
library(here)

# ---- 1. Configuration -----------------------------------------------------------

models <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "HCLIM_ERA5",
            "RACMO_CESM2", "RACMO_MPI_ESM1", "RACMO_ERA5",
            "MetUM_ERA5")

input_base  <- here("Data/Environmental_predictors/PolarRes26")
output_base <- here("Data/Environmental_predictors/PolarRes26/Regridded")

# Target-grid template: extent/CRS/resolution ONLY -- never used to mask.
ant_extent_grid_1km_path <- here("Data/Ant_extent_grid_1km.tif")

# Ice-free mask, applied to every output regardless of historical/future.
ice_free_domain_1km_path <- here("Data/ice_free_domain_1km.tif")

ant_extent_grid_1km <- rast(ant_extent_grid_1km_path)
ice_free_domain_1km <- rast(ice_free_domain_1km_path)

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

# If TRUE, skip an input file entirely (no reprojection/resampling work at
# all) when every output it would produce already exists on disk. Set to
# FALSE to always reprocess and overwrite.
SKIP_EXISTING <- TRUE

# ---- 2. Build the job table: every .tif in each model's "comparison" folder ----

job_table <- map_dfr(models, function(model) {
  model_dir <- file.path(input_base, model, "comparison")
  if (!dir.exists(model_dir)) {
    message("Comparison directory not found, skipping: ", model_dir)
    return(tibble())
  }
  files <- list.files(model_dir, pattern = "\\.tif$", full.names = TRUE)
  tibble(model = model, input_path = files, filename = basename(files))
})

if (nrow(job_table) == 0) {
  stop("No .tif files found under any <model>/comparison folder under ", input_base,
       " -- has the comparison version of Script 2 been run yet?")
}

message(nrow(job_table), " output file(s) found across ",
        length(unique(job_table$model)), " model(s).")

# ---- 3. job_index selects ONE input file -----------------------------------------

args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])

this_job   <- job_table[job_index, ]
model      <- this_job$model
input_path <- this_job$input_path
filename   <- this_job$filename
message("Job ", job_index, "/", nrow(job_table), " -> ", model, " / ", filename)

# ---- 4. Output paths -------------------------------------------------------------
# Same "comparison" subfolder, but under Regridded/<model>/ instead of
# <model>/ directly.

out_dir <- file.path(output_base, model, "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stem <- tools::file_path_sans_ext(filename)
icefree_path <- file.path(out_dir, paste0(stem, "_ICEFREE.tif"))

expected_outputs <- c(icefree_path)

# ---- 4b. Skip-existing check ------------------------------------------------------

if (SKIP_EXISTING && all(file.exists(expected_outputs))) {
  message("  all expected output(s) already exist -- skipping: ",
          paste(basename(expected_outputs), collapse = ", "))
  quit(save = "no", status = 0)
}

# ---- 5. Project, resample, mask, save --------------------------------------------

r <- rast(input_path)

# Hardcode the resolution to 11km
reso <- 11000

if (is.na(crs(r)) || crs(r) == "") {
  
  # Known HCLIM quirk: some HCLIM outputs come out of Script 1/2 with no
  # CRS attached, even though their grid/extent is fine. For HCLIM only,
  # recover the CRS from a known-good HCLIM template file -- but only after
  # confirming the extents actually line up, so we never silently stamp a
  # wrong CRS onto a raster that doesn't actually match the template grid.
  if (grepl("^HCLIM", model)) {
    
    hclim_template <- rast(hclim_crs_template_path)
    
    extents_match <- isTRUE(all.equal(
      as.vector(ext(r)), as.vector(ext(hclim_template)),
      tolerance = 1e-6
    ))
    
    if (!extents_match) {
      stop("Input raster has no CRS and its extent does not match the HCLIM ",
           "CRS template (", hclim_crs_template_path, ") -- refusing to guess ",
           "a CRS for a grid that doesn't line up. Input: ", input_path,
           "\n  input extent:    ", paste(round(as.vector(ext(r)), 4), collapse = ", "),
           "\n  template extent: ", paste(round(as.vector(ext(hclim_template)), 4), collapse = ", "))
    }
    
    crs(r) <- crs(hclim_template)
    message("  input had no CRS -- extent matched the HCLIM template, so CRS was ",
            "copied from it: ", input_path)
    
  } else {
    # RACMO / MetUM (or anything else): an empty CRS here is not a known,
    # safe-to-patch quirk -- treat it as the hard error it always has been.
    stop("Input raster has no CRS: ", input_path,
         " -- reprojecting without a source CRS would silently produce garbage.",
         " Check this file (and Script 1/2's terra version, if this is a RACMO file).")
  }
}

# Step 1: Make a domain for target CRS with matching resolution
target_grid <- rast(extent = ext(ant_extent_grid_1km), crs = crs(ant_extent_grid_1km),
                    resolution = reso)

# Sanity check: extent must be a clean integer multiple of the resolution,
# otherwise project()/rast() will silently produce a grid that doesn't
# tile evenly and downstream resample() will misalign.
stopifnot(
  (xmax(target_grid) - xmin(target_grid)) %% res(target_grid) == 0,
  (ymax(target_grid) - ymin(target_grid)) %% res(target_grid) == 0
)

# Step 2: Reproject into the domain's CRS with bilinear interpolation
r <- project(r, target_grid, method = "bilinear")

# Step 3: Resample to the 1km domain grid (already in the right CRS &
# extent, so "near" here is just chopping up 11km to 1km grids)
r <- resample(r, ant_extent_grid_1km, method = "near")

# ---- 6. Mask, save ----------------------------------------------------------------
# No coastline cropping and no sea-ice special case here -- every file just
# gets the single ice-free-masked output.

r_icefree <- mask(r, ice_free_domain_1km, maskvalue = NA)
writeRaster(r_icefree, icefree_path, overwrite = TRUE)
message("  wrote ", icefree_path)

message("Done: ", model, " / ", filename)