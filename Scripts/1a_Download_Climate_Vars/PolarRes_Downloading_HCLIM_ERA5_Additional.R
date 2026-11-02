# HPC version

library(terra)
library(glue)
library(gdalraster)
library(vapour)

# HCLIM_ERA5 -- ADDITIONAL VARIABLES (rsus, siconca, snm)
# --------------------------------------------------------------------------
#   - no reprojection / no .tif conversion -- the .nc files are kept as-is
#   - everything goes into one flat "additional" folder
#

u <- "url" # This was a folder shared with the authors by the PolarRES team - unable to be supplied here

model_name <- "HCLIM_ERA5_"

dirpath <- "/Data/PolarRes26"
outdir  <- glue("{dirpath}/{model_name}")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Listing of the ERA5 folder (netCDF files only)
files  <- vsi_read_dir(glue("/vsicurl/{u}"))
files  <- files[grepl("\\.nc$", files)]
df_all <- data.frame(file = files, url = glue("{u}{files}"))

# Variables to download:
#   rsus    - surface upwelling shortwave radiation
#   siconca - sea ice concentration (%)
#   snm     - snowmelt (kg m-2 s-1)
vars <- c("rsus", "siconca", "snm")

for (v in vars) {
  # match on the filename prefix (e.g. "rsus_fp_ANT11_...") so a variable
  # name appearing elsewhere in a filename can't cause a false match
  df <- df_all[grepl(glue("^{v}_"), df_all$file), ]
  message(glue("{v}: {nrow(df)} files found"))
  
  for (i in seq_len(nrow(df))) {
    localfile <- glue("{outdir}/{df$file[i]}")
    if (!file.exists(localfile)) {
      # download to a temp name, then rename, so an interrupted job
      # doesn't leave a truncated .nc that gets skipped on re-run
      tmpfile <- paste0(localfile, ".part")
      ok <- tryCatch({
        curl::curl_download(df$url[i], tmpfile)
        file.rename(tmpfile, localfile)
        TRUE
      }, error = function(e) {
        message(glue("FAILED: {df$file[i]} -- {conditionMessage(e)}"))
        if (file.exists(tmpfile)) file.remove(tmpfile)
        FALSE
      })
    }
  }
}