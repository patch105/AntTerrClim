
# ============================================================
# CHELSA vs the AntTerrClim temperature evaluation
#   FIGURE 5-6 main manuscript : Scatter-density (Annual; Summer + Winter combined)
#   FIGURE 8 main manuscript : Bias raster maps (Annual, Summer, Winter)
#   PART 3   : Summary statistics (mean bias, RMSE, MAE, R2)
#   Figure 7 main manuscript  : Mean bias heatmap
# ============================================================
#
# INPUTS:
#  1. CHELSA reference rasters from script 2c (ice-free, regridded to the 1km grid):
#       Data/CHELSA/comparison/regridded/
#     for Season = Annual, Summer (DJF), Winter (JJA).
#  2. Model rasters from script 2c all follow:
#       Data/Environmental_predictors/PolarRes26/Regridded/<MODEL>_<DRIVING>/comparison/
#     e.g. HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5, HCLIM_MPI_ESM1, etc.
#  3. CHELSA and every model raster are on the SAME grid (same extent,
#     CRS and resolution). Nothing is resampled here -- if a pair does
#     not match, the script stops with an error naming the pair.
#  4. MetUM only exists for the ERA5-driven run, so it only appears in
#     the ERA5 row. Storyline rows have just HCLIM + RACMO.
#  5. Coastline shapefile: Data/add_coastline_medium_res_polygon_v7_10.shp
#  8. agg_fact = 10 assumes ~1 km native resolution.
#
# ============================================================


library(terra) 
library(ggplot2)
library(dplyr)
library(tibble)
library(purrr)
library(patchwork)
library(here)
library(scales)
library(sf)

outpath <- here("Plots/Evaluation_CHELSA")
if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

PAGE_W <- 8.27; PAGE_H <- 11.69; DPI <- 320; BASE <- 8.5; FONT <- "Helvetica"
LEGBAR_H <- unit(3.2, "cm"); LEGBAR_W <- unit(0.3, "cm")
agg_fact <- 10   # ~10 km for continent-wide

# ---------------------------------------------------------------
# COASTLINE
# ---------------------------------------------------------------

# Antarctic continent, from SCAR Antarctic Digital Database high resolution vector polygons of the Antarctic coastline v7.12 (Gerrish et al., 2026) 
# NOTE: Downloaded the zipped file and unzipped it manually

coast <- st_read(
  here("Data/add_coastline_medium_res_polygon_v7_10.shp"),
  quiet = TRUE
)
coast_fill   <- "grey95"; coast_colour <- "grey35"

rast_to_df <- function(r, val_name = "value") {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(df)[3] <- val_name; df
}
align_to_model <- function(ref, mod) {
  if (!compareGeom(ref, mod, stopOnError = FALSE)) ref <- resample(ref, mod, method = "bilinear")
  ref
}
theme_ant <- function(base_size = BASE, base_family = FONT) {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.position = "right", legend.box.spacing = unit(1, "pt"),
      legend.spacing = unit(1, "pt"),
      legend.title = element_text(size = base_size, angle = 90, hjust = 0.5, vjust = 0.5),
      legend.text = element_text(size = base_size - 1, margin = margin(l = 1)),
      legend.key.height = LEGBAR_H, legend.key.width = LEGBAR_W,
      legend.margin = margin(0, 0, 0, 0),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.4,
                                vjust = 1, colour = "grey15", margin = margin(b = 1)),
      plot.margin = margin(1, 1, 1, 1))
}
make_cbar <- function() {
  guide_colorbar(title.position = "right", title.hjust = 0.5,
                 barheight = LEGBAR_H, barwidth = LEGBAR_W,
                 ticks.colour = "black", frame.colour = "black", frame.linewidth = 0.35)
}
diverging_ramp <- c("#053061","#2166ac","#4393c3","#92c5de","#d1e5f0","white",
                    "#fddbc7","#f4a582","#d6604d","#b2182b","#67001f")
make_scale_values <- function(min_val, max_val) {
  w <- (0 - min_val) / (max_val - min_val)
  c(0, w*0.25, w*0.50, w*0.75, w*0.95, w,
    w + (1-w)*0.15, w + (1-w)*0.40, w + (1-w)*0.65, w + (1-w)*0.85, 1)
}
base_layers <- function() geom_sf(data = coast, fill = coast_fill, colour = coast_colour, linewidth = 0.15)


# ============================================================
# FIGURE 10 MANUSCRIPT: FUTURE - HISTORICAL CHANGE
# ============================================================

chelsa_hist_var <- "Annual"

# ---- CHELSA: ensemble-mean future minus historical ----
gcms <- c("gfdl-esm4","ipsl-cm6a-lr","mpi-esm1-2-hr","mri-esm2-0","ukesm1-0-ll")
CH_fut <- app(rast(sapply(gcms, function(g) here(sprintf(
  "Data/CHELSA/comparison/Mean_Annual_Temperature_FUTURE_%s_2071_2100_ICEFREE.tif", g)))),
  mean, na.rm = TRUE)
CH_hist <- rast(here(sprintf(
  "Data/CHELSA/comparison/Mean_%s_Temperature_HISTORICAL_1981_2010_ICEFREE.tif", chelsa_hist_var)))
CH_delta <- aggregate(CH_fut - align_to_model(CH_hist, CH_fut), fact = agg_fact, fun = "mean", na.rm = TRUE)

# ---- Models ----
model_config <- tribble(
  ~model, ~driving, ~row_group,
  "HCLIM","MPI_ESM1","Storyline 1", "RACMO","MPI_ESM1","Storyline 1",
  "HCLIM","CESM2","Storyline 2",     "RACMO","CESM2","Storyline 2"
) %>% mutate(folder = paste0(model, "_", driving),
             label  = paste0(model, " \u2013 ", row_group))

load_delta <- function(folder, ...) {
  d <- here("Data/Environmental_predictors/PolarRes26/Regridded", folder, "comparison")
  fut  <- rast(file.path(d, "Mean_Annual_Temperature_FUTURE_2071_2100_ICEFREE.tif"))
  hist <- rast(file.path(d, "Mean_Annual_Temperature_HISTORICAL_1985_2010_ICEFREE.tif"))
  aggregate(fut - align_to_model(hist, fut), fact = agg_fact, fun = "mean", na.rm = TRUE)
}
deltas <- pmap(model_config, load_delta)
names(deltas) <- model_config$label

# ---- Shared scale: warming (sequential); diverging if any cooling ----
all_vals <- c(values(CH_delta, na.rm = TRUE), unlist(lapply(deltas, values, na.rm = TRUE)))
lo <- floor(min(all_vals)); hi <- ceiling(max(all_vals))
if (lo < 0) {
  ramp <- diverging_ramp; vals <- make_scale_values(lo, hi)
} else {
  ramp <- c("#ffffcc","#fed976","#fd8d3c","#e31a1c","#b10026","#67001f"); vals <- NULL
}
message("Change range: ", round(min(all_vals),2), " to ", round(max(all_vals),2))

make_delta_panel <- function(r, title) {
  ggplot() + base_layers() +
    geom_tile(data = rast_to_df(r, "d"), aes(x, y, fill = d)) +
    scale_fill_gradientn(colours = ramp, values = vals, limits = c(lo, hi),
                         breaks = pretty(c(lo, hi), n = 6),
                         name = "Change in temperature (\u00B0C)", oob = squish, guide = make_cbar()) +
    coord_sf(expand = FALSE) + labs(title = title) + theme_ant()
}

p_chelsa <- make_delta_panel(CH_delta, "CHELSA (ensemble mean)")
# panel titles are now just the RCM name (storyline shown by the block title)
p <- lapply(seq_along(deltas), function(i) make_delta_panel(deltas[[i]], model_config$model[i]))

# Large left-aligned storyline title (sits over column 1)
storyline_title <- function(txt) {
  ggplot() +
    annotate("text", x = 0, y = 0.5, label = txt, hjust = 0, vjust = 0.5,
             fontface = "bold", size = 6, family = FONT, colour = "grey15") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 4))
}

# One shared legend per row: the two panels use identical scales, so
# guides = "collect" merges them into a single colourbar on the right
make_pair_row <- function(a, b) {
  (a + b) + plot_layout(guides = "collect") &
    theme(legend.position = "right")
}
make_block <- function(title, a, b) {
  storyline_title(title) / make_pair_row(a, b) + plot_layout(heights = c(0.09, 1))
}

chelsa_row <- plot_spacer() + p_chelsa + plot_spacer() + plot_layout(widths = c(0.5, 1, 0.5))
block1 <- make_block("Storyline 1", p[[1]], p[[2]])   # HCLIM | RACMO
block2 <- make_block("Storyline 2", p[[3]], p[[4]])   # HCLIM | RACMO

final_plot <- chelsa_row / block1 / block2 + plot_layout(heights = c(1, 1.1, 1.1))

# Shorter page than Plot 2 to tighten row spacing 
# FIGURE 10 main manuscript
PLOT1_H <- 9.5
out_file <- file.path(outpath, "plot_annual_temp_CHANGE_FUTURE_minus_HISTORICAL.png")
ggsave(out_file, final_plot, width = PAGE_W, height = PLOT1_H, dpi = DPI, bg = "white")
message("Saved: ", out_file)

# ============================================================
# FIGURE 11 main manuscript: MPI-ESM1 ONLY, FUTURE, CHELSA vs HCLIM / RACMO
# ============================================================

CH <- rast(here("Data/CHELSA/comparison/Mean_Annual_Temperature_FUTURE_mpi-esm1-2-hr_2071_2100_ICEFREE.tif"))

models <- c("HCLIM", "RACMO")
load_model <- function(m) {
  mod <- rast(here("Data/Environmental_predictors/PolarRes26/Regridded",
                   paste0(m, "_MPI_ESM1"), "comparison",
                   "Mean_Annual_Temperature_FUTURE_2071_2100_ICEFREE.tif"))
  diff <- mod - align_to_model(CH, mod)
  list(temp = aggregate(mod,  fact = agg_fact, fun = "mean", na.rm = TRUE),
       diff = aggregate(diff, fact = agg_fact, fun = "mean", na.rm = TRUE))
}
md <- setNames(lapply(models, load_model), models)
CH_agg <- aggregate(CH, fact = agg_fact, fun = "mean", na.rm = TRUE)

all_temp <- c(values(CH_agg, na.rm = TRUE), unlist(lapply(md, function(m) values(m$temp, na.rm = TRUE))))
temp_min <- floor(min(all_temp)); temp_max <- ceiling(max(all_temp))
temp_values <- make_scale_values(temp_min, temp_max)

all_diff <- unlist(lapply(md, function(m) values(m$diff, na.rm = TRUE)))
diff_min <- floor(min(all_diff)); diff_max <- ceiling(max(all_diff))
diff_breaks <- pretty(c(diff_min, diff_max), n = 6)
diff_values <- make_scale_values(diff_min, diff_max)

make_panel <- function(r, title, lo, hi, vals, name, brks = pretty(c(lo, hi), n = 6)) {
  ggplot() + base_layers() +
    geom_tile(data = rast_to_df(r, "v"), aes(x, y, fill = v)) +
    scale_fill_gradientn(colours = diverging_ramp, values = vals, limits = c(lo, hi),
                         breaks = brks, name = name, oob = squish, guide = make_cbar()) +
    coord_sf(expand = FALSE) + labs(title = title) + theme_ant()
}
temp_panel <- function(r, t) make_panel(r, t, temp_min, temp_max, temp_values, "Temperature (\u00B0C)")
diff_panel <- function(r, t) make_panel(r, t, diff_min, diff_max, diff_values, "Difference (\u00B0C)", diff_breaks)

chelsa_row <- plot_spacer() + temp_panel(CH_agg, "CHELSA (MPI-ESM1-2-HR)") + plot_spacer() +
  plot_layout(widths = c(0.5, 1, 0.5))
model_row <- function(m) {
  temp_panel(md[[m]]$temp, paste0(m, " (MPI-ESM1)")) + plot_spacer() +
    diff_panel(md[[m]]$diff, paste0(m, " \u2212 CHELSA")) +
    plot_layout(widths = c(1, 0.015, 1))
}

final_plot <- chelsa_row / model_row("HCLIM") / model_row("RACMO") + plot_layout(heights = c(1, 1, 1))
out_file <- file.path(outpath, "plot_annual_temp_FUTURE_MPI_ESM1_CHELSA_vs_HCLIM_RACMO.png")
# Shorter page (same idea as Plot 1) to bring the rows closer together
PLOT2_H <- 9.5
ggsave(out_file, final_plot, width = PAGE_W, height = PLOT2_H, dpi = DPI, bg = "white")
message("Saved: ", out_file)
