# HPC version

library(terra)
library(glue)
library(gdalraster)
library(vapour)

# HCLIM_CESM2 -- ADDITIONAL VARIABLES (rsus, snm, snc, siconca)
# --------------------------------------------------------------------------
# Same source/style as the original HCLIM_MPI_ESM1 download script, but:
#   - no reprojection / no .tif conversion -- the .nc files are kept as-is
#   - everything goes into one flat "additional" folder, not split into
#     historical/ssp370 subfolders (that can be sorted out later -- for now
#     this just gets every matching file downloaded)
#

u <- "url" # This was a folder shared with the authors by the PolarRES team - unable to be supplied here

model_name <- "HCLIM_CESM2_"

dirpath <- "/Data/PolarRes26"
outdir  <- glue("{dirpath}/{model_name}additional")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Combined listing across all three source folders, restricted to files
# that are actually CESM2's (MPI's files live in some of the same folders)
df_list <- list()
for (u_i in u) {
  files <- vsi_read_dir(glue("/vsicurl/{u_i}"))
  df_list[[u_i]] <- data.frame(file = files, url = glue("{u_i}{files}"))
}
df_all <- do.call(rbind, df_list)
df_all <- df_all[grepl("CESM2", df_all$file), ]

# NOTE: since CESM2 is split across three folders, the same filename could
# in principle exist in more than one of them (the subset diagnostics
# script already found exactly this for one snm file, split between CESM/
# and midfut/) -- file.exists() below skips anything already downloaded, so
# only the first copy encountered (in the order `u` is listed above) is
# kept, and a later duplicate is left alone rather than re-downloaded or
# overwritten. Worth re-running the subset diagnostics on this "additional"
# folder afterwards to check which filename(s) that affected.


# RSUS --------------------------------------------------------------------

df <- df_all[grepl("rsus", df_all$file), ]

for (i in seq_len(nrow(df))) {
  localfile <- glue("{outdir}/{df$file[i]}")
  if (!file.exists(localfile)) {
    curl::curl_download(df$url[i], localfile)
  }
}


# SNM (kg m-2 s-1) --------------------------------------------------------------------

df <- df_all[grepl("snm", df_all$file), ]

for (i in seq_len(nrow(df))) {
  localfile <- glue("{outdir}/{df$file[i]}")
  if (!file.exists(localfile)) {
    curl::curl_download(df$url[i], localfile)
  }
}


# SNC (% "Snow area percentage") --------------------------------------------------------------------

df <- df_all[grepl("snc", df_all$file), ]

for (i in seq_len(nrow(df))) {
  localfile <- glue("{outdir}/{df$file[i]}")
  if (!file.exists(localfile)) {
    curl::curl_download(df$url[i], localfile)
  }
}


# SICONCA (%) --------------------------------------------------------------------

df <- df_all[grepl("siconca", df_all$file), ]

for (i in seq_len(nrow(df))) {
  localfile <- glue("{outdir}/{df$file[i]}")
  if (!file.exists(localfile)) {
    curl::curl_download(df$url[i], localfile)
  }
}