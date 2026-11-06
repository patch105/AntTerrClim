# ==============================================================================
# DOWNLOAD THE RAW CLIMATE VARIABLES
# ==============================================================================

# This script will download the climate data from the PolarRES project for a subset of variables and timeframes

# Designed to be run on a High-Performance Computer with an array of jobs

# The script creates a folder directory from the 'outbase' location that matches the directory on the PolarRES server, saving each .nc file within a given directory

# Since the files to download are quite large, this script takes a couple of hours to run


library(arrow)
library(dplyr)
library(stringr)
library(glue)

model_name <- "MetUM_ERA5"
model_full <- "ERA5"

# Output base directory
outbase <- "/Data/PolarRes26"

# ── Load parquet ─────────────────────────────────────────────────────────────
data_source <- arrow::read_parquet(
  "https://github.com/mdsumner/aad-filelist/releases/download/latest/cordex-polarres-ant12-metum-bas.parquet"
)


# ── Variables to check ───────────────────────────────────────────────────────
variables <- c( # Removed snm for now since it's incorrect
  "tas",
  "tasmax",
  "tasmin",
  "sfcWind",
  "siconca",
  "pr",
  "snc",
  "rsus",
  "rsds",
  "hurs"
)

# ── Job index from command line ───────────────────────────────────────────────
args      <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])
variable  <- variables[[job_index]]
print(paste0("Variable is: ", variable))

# ── Subset URLs ───────────────────────────────────────────────────────────────
# Keep only:
#   1. Files containing the target variable in their path (e.g. .../day/tas/...)
#   2. Daily frequency only (path segment "day" immediately before the variable)
#   3. ERA5-driven run (ERA5 appears as a token in the filename, e.g.
#      tas_ANT-12_ERA5_evaluation_r1i1p1f1_BAS_MetUM_v1-r1_day_....nc)
dat <- data_source %>%
  filter(
    # variable folder in URL path
    str_detect(url, paste0("/day/", variable, "/")),
    str_detect(url, model_full)
  )

if (nrow(dat) == 0) {
  message(sprintf("No daily files found for variable '%s' / ERA5. Exiting.", variable))
  quit(status = 0)
}

message(sprintf("Found %d file(s) to download for %s.", nrow(dat), variable))

# ── Build local paths ─────────────────────────────────────────────────────────
# New URL structure (JASMIN GWS access), e.g.:
#   https://gws-access.jasmin.ac.uk/public/polarres/CORDEX-output/v20260129/
#     ANT-12/day/tas/tas_ANT-12_ERA5_evaluation_r1i1p1f1_BAS_MetUM_v1-r1_day_19950101-19951231.nc
#
# Unlike the old THREDDS layout, there's no "evaluation/" (or "sspXXX/") *folder*
# segment to key off -- the scenario/experiment info is embedded in the filename
# itself. The directory structure that remains is just ANT-12/<freq>/<var>/<file>.
#
# We want local path:
#   <outbase>/<model_name>/ANT-12/<freq>/<var>/<file>.nc
#
# Strategy: capture everything from "ANT-12/" onward in the URL and place it
# under <outbase>/<model_name>/
dat <- dat %>%
  mutate(
    # Capture the relative path from "ANT-12/" to the end
    rel_path = str_extract(url, "ANT-12/.+$"),
    local_path = file.path(outbase, model_name, rel_path)
  )

# ── Download loop ─────────────────────────────────────────────────────────────
for (i in seq_len(nrow(dat))) {
  
  url_i   <- dat$url[i]
  local_i <- dat$local_path[i]
  
  # Create directory if needed
  dir.create(dirname(local_i), recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(local_i)) {
    message(sprintf("[%d/%d] Already exists, skipping: %s", i, nrow(dat), basename(local_i)))
    next
  }
  
  message(sprintf("[%d/%d] Downloading: %s", i, nrow(dat), basename(local_i)))
  
  tryCatch(
    curl::curl_download(url_i, local_i),
    error = function(e) {
      message(sprintf("  ERROR downloading %s: %s", basename(url_i), conditionMessage(e)))
      # Remove partial file if download failed
      if (file.exists(local_i)) file.remove(local_i)
    }
  )
}

message("Done.")