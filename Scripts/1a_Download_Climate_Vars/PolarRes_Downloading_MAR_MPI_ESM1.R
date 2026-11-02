# HPC version

library(arrow)
library(dplyr)
library(stringr)
library(glue)

model_name <- "MAR_MPI_ESM1"
model_full <- "MPI-ESM1-2-HR"

# Output base directory 
outbase <- "/Data/PolarRes26"

# ── Load parquet ─────────────────────────────────────────────────────────────
data_source <- arrow::read_parquet(
  "https://github.com/mdsumner/aad-filelist/releases/download/latest/polarres-ant12-mar-v3.13-ulg.parquet"
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
#   3. MPI-ESM1-2-LR GCM (MPI_ESM1 in model_name maps to MPI-ESM1-2-LR in URLs)

dat <- data_source %>%
  filter(
    # variable folder in URL path
    str_detect(url, paste0("/day/", variable, "/")),
    str_detect(url, model_full)
  )

if (nrow(dat) == 0) {
  message(sprintf("No daily files found for variable '%s' / MPI-ESM1-2-LR. Exiting.", variable))
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