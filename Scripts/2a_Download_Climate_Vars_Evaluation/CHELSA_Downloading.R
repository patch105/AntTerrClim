
# ==============================================================================
# Download the FUTURE CHELSA variables for comparison with AntTerrClim 
# ==============================================================================

# HPC version - CHELSA temperature climatologies
# Historical data for 1981-2010
# SSP370, 2071-2100, models: GFDL-ESM4, IPSL-CM6A-LR, MRI-ESM2-0, MPI-ESM1-2-HR, UKSEM1-0-LL
# Monthly (01-12) GeoTIFF files 

# Designed to be run on a High-Performance Computer with an array of jobs

# Every job is a separate input file 

# The script creates a folder directory structured as: here("Data/CHELSA") where the CHELSA data are saved as .tif files

library(terra)
library(glue)
library(gdalraster)
library(vapour)
library(here)

# Extract arguments from command line
args <- commandArgs(trailingOnly = TRUE)
job_index <- as.integer(args[1])

# Stagger jobs to avoid hammering the server simultaneously
Sys.sleep(runif(1, min = job_index * 2, max = job_index * 2 + 30))

# Full URL list (60 files)
base <- "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/climatologies/tas"

months <- sprintf("%02d", 1:12)

future_models <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
future_slugs  <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0",  "ukesm1-0-ll")

# Build future URLs: model x month (model varies slowest)
future_urls <- unlist(lapply(seq_along(future_models), function(m) {
  glue("{base}/2071-2100/{future_models[m]}/ssp370/",
       "CHELSA_{future_slugs[m]}_r1i1p1f1_w5e5_ssp370_tas_{months}_2071-2100_V.2.1.tif")
}))

# Build historical URLs
hist_urls <- glue("{base}/1981-2010/CHELSA_tas_{months}_1981-2010_V.2.1.tif")

# All 60 URLs in order (jobs 1-48 = future, 49-60 = historical)
all_urls <- c(future_urls, hist_urls)

# Select file for this job (1-based)
file_url  <- all_urls[job_index]
file_name <- basename(file_url)

dirpath   <- here("Data/CHELSA")
localfile <- file.path(dirpath, file_name)

# --- Check: file already exists ---
if (file.exists(localfile)) {
  message(glue("Job {job_index}: Already exists — skipping {file_name}"))
  quit(status = 0)
}

message(glue("Job {job_index}: Downloading {file_name}"))

max_attempts <- 5
for (attempt in 1:max_attempts) {
  result <- tryCatch({
    curl::curl_download(file_url, localfile)
    "success"
  }, error = function(e) {
    message(glue("Attempt {attempt} failed: {e$message}"))
    "error"
  })
  
  if (result == "success") {
    message(glue("Download complete: {localfile}"))
    break
  } else if (attempt < max_attempts) {
    wait <- 30 * attempt
    message(glue("Retrying in {wait}s..."))
    Sys.sleep(wait)
  } else {
    stop(glue("Failed to download {file_name} after {max_attempts} attempts"))
  }
}
