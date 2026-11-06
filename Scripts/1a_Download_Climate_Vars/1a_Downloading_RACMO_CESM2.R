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

model_name <- "RACMO_CESM2"
model_full <- "CESM2"

# Output base directory 
outbase <- "/Data/PolarRes26"

# ── Load parquet ─────────────────────────────────────────────────────────────
data_source <- arrow::read_parquet(
  "https://github.com/mdsumner/aad-filelist/releases/download/latest/cordex-polarres-ant12-racmo-uu-imau.parquet"
)

# ── Variables to check ───────────────────────────────────────────────────────
variables <- c(
  "tas",
  "tasmax",
  "tasmin",
  "sfcWind",
  "siconca",
  "pr",
  "snc",
  "snm",
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

dat <- data_source %>%
  filter(
    # variable folder in URL path
    str_detect(url, paste0("/day/", variable, "/")),
    str_detect(url, model_full)
  )

if (nrow(dat) == 0) {
  message(sprintf("No daily files found for variable '%s' / . Exiting.", variable))
  quit(status = 0)
}

message(sprintf("Found %d file(s) to download for %s.", nrow(dat), variable))

# ── Build local paths ─────────────────────────────────────────────────────────
# URL structure (after the THREDDS prefix):
#   .../PolarRes/ANT-12/UU-IMAU/<GCM>/<scenario>/.../<freq>/<var>/.../<file>.nc
#
# We want local path:
#   <outbase>/<model_name>/<scenario>/.../<freq>/<var>/.../<file>.nc
#
# Strategy: extract everything from "historical" or "ssp..." onward in the URL
# and place it under <outbase>/<model_name>/

dat <- dat %>%
  mutate(
    # Capture the relative path from "historical" or "sspXXX" to the end
    rel_path = str_extract(url, "((?:historical|ssp[^/]+)/.+)"),
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