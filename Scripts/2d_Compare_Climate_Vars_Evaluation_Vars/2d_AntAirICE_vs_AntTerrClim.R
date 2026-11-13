# ============================================================
# AntAirICE vs the AntTerrClim temperature evaluation
#   FIGURE 5-6 main manuscript : Scatter-density (Annual; Summer + Winter combined)
#   FIGURE 8 main manuscript : Bias raster maps (Annual, Summer, Winter)
#   PART 3   : Summary statistics (mean bias, RMSE, MAE, R2)
#   Figure 7 main manuscript  : Mean bias heatmap
# ============================================================
#
# INPUTS:
#  1. AntAirICE reference rasters from script 2c (ice-free, regridded to the 1km grid):
#       Data/AntAirICE/comparison/regridded/
#         Mean_<Season>_Temperature_HISTORICAL_2003_2014_ICEFREE.tif
#     for Season = Annual, Summer (DJF), Winter (JJA).
#  2. Model rasters from script 2c all follow:
#       Data/Environmental_predictors/PolarRes26/Regridded/<MODEL>_<DRIVING>/comparison/
#         Mean_<Season>_Temperature_HISTORICAL_2003_2014_ICEFREE.tif
#     e.g. HCLIM_ERA5, RACMO_ERA5, MetUM_ERA5, HCLIM_MPI_ESM1, etc.
#  3. AntAirICE and every model raster are on the SAME grid (same extent,
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
library(tidyr)
library(purrr)
library(patchwork)
library(here)
library(viridis)
library(Metrics)
library(sf)
library(scales)

outpath <- here("Plots/Evaluation_AntAirICE")
if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

# ---------------------------------------------------------------
# SHARED CONFIGURATION
# ---------------------------------------------------------------
seasons <- c("Annual", "Summer", "Winter")

ref_path <- function(season) {
  here(
    "Data/AntAirICE/comparison/regridded",
    paste0("Mean_", season, "_Temperature_HISTORICAL_2003_2014_ICEFREE.tif")
  )
}

model_path <- function(folder, season) {
  here(
    "Data/Environmental_predictors/PolarRes26/Regridded",
    folder, "comparison",
    paste0("Mean_", season, "_Temperature_HISTORICAL_2003_2014_ICEFREE.tif")
  )
}

# AntAirICE and the models are expected to share a grid. Stop loudly
# (rather than silently resampling) if they don't.
check_same_grid <- function(ref, mod, label = "") {
  if (!isTRUE(compareGeom(ref, mod, stopOnError = FALSE))) {
    stop("AntAirICE and model raster are not on the same grid",
         if (nzchar(label)) paste0(" (", label, ")") else "",
         " -- check extent / CRS / resolution of the regridded files.")
  }
  invisible(TRUE)
}

extract_pairs <- function(ref, mod, label) {
  check_same_grid(ref, mod, label)
  vals <- data.frame(
    x = as.vector(values(ref)),
    y = as.vector(values(mod))
  )
  vals <- vals[complete.cases(vals), ]
  vals$model <- label
  vals
}


# ###############################################################
# FIGURE 1: SCATTER-DENSITY (Storyline 1 / Storyline 2 / ERA5)
# ###############################################################

# ---------------------------------------------------------------
# STYLE
# ---------------------------------------------------------------
row_title_colour <- "grey35"
hide_inner_y_labels <- TRUE # only first panel in a row shows y tick labels
row_height_pad   <- 0.10    # extra relative height per row for title / x-label
# (increase if rows look cramped, decrease if gaps)

# Shared axis titles (one x + one y per row)
x_axis_title <- "Temperature (°C) AntAir ICE"
y_axis_title <- "Temperature (°C) Model"

theme_storyline <- function(base_size = 11, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background = element_blank(),
      # margin(b = ...) lifts the panel titles (HCLIM etc.) off the plot box
      strip.text       = element_text(face = "bold", size = 10.5, colour = "grey20",
                                      margin = margin(b = 5)),
      axis.title       = element_text(size = 9.5, colour = "grey30"),
      axis.text        = element_text(size = 8.5, colour = "grey45"),
      axis.line        = element_line(colour = "grey60", linewidth = 0.35),
      axis.ticks       = element_line(colour = "grey60", linewidth = 0.35),
      panel.border     = element_rect(colour = "grey60", fill = NA, linewidth = 0.35),
      panel.spacing    = unit(0.45, "cm"),
      legend.title     = element_text(size = 9.5, colour = "grey20"),
      legend.text      = element_text(size = 8.5, colour = "grey40")
    )
}

# A plain text "plot" used for row titles and shared axis titles
label_plot <- function(label, angle = 0, size = 3.4, face = "plain",
                       colour = "grey30", x = 0.5, hjust = 0.5) {
  ggplot() +
    annotate("text", x = x, y = 0.5, label = label, angle = angle,
             size = size, fontface = face, colour = colour,
             hjust = hjust, family = "Helvetica") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void()
}

# ---------------------------------------------------------------
# MODEL CONFIGURATION
# ---------------------------------------------------------------
row_levels <- c("Storyline 1", "Storyline 2", "ERA5 Re-analysis")

model_config <- tribble(
  ~model,   ~driving,     ~row_group,
  "HCLIM",  "MPI_ESM1",   "Storyline 1",
  "RACMO",  "MPI_ESM1",   "Storyline 1",
  "HCLIM",  "CESM2",      "Storyline 2",
  "RACMO",  "CESM2",      "Storyline 2",
  "HCLIM",  "ERA5",       "ERA5 Re-analysis",
  "RACMO",  "ERA5",       "ERA5 Re-analysis",
  "MetUM",  "ERA5",       "ERA5 Re-analysis"
) %>%
  mutate(
    col_label = model,
    folder    = paste0(model, "_", driving)
  )

# ---------------------------------------------------------------
# PAIRED DATA
# ---------------------------------------------------------------
build_season_pairs <- function(season) {
  ref <- rast(ref_path(season))
  
  pmap_dfr(model_config, function(model, driving, row_group, col_label, folder) {
    mod <- rast(model_path(folder, season))
    df  <- extract_pairs(ref, mod, label = paste(model, driving, sep = "_"))
    df$row_group <- row_group
    df$col_label <- col_label
    df$season    <- season
    df
  })
}

get_shared_fill_limit <- function(df, breaks) {
  x_bin <- cut(df$x, breaks = breaks, include.lowest = TRUE)
  y_bin <- cut(df$y, breaks = breaks, include.lowest = TRUE)
  grp   <- interaction(df$row_group, df$col_label, df$season, drop = TRUE)
  max(table(grp, x_bin, y_bin))
}

# ---------------------------------------------------------------
# PANELS (no axis titles - those are shared per row)
# ---------------------------------------------------------------
build_panel_plot <- function(df_panel, strip_label, ax_lim, ax_breaks,
                             common_breaks, fill_limit, show_y_text = TRUE) {
  df_panel$strip <- strip_label
  
  p <- ggplot(df_panel, aes(x = x, y = y)) +
    stat_bin2d(breaks = list(x = common_breaks, y = common_breaks),
               aes(fill = after_stat(count))) +
    geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.4) +
    scale_fill_viridis_c(
      option   = "plasma",
      name     = "Number of\ngrid cells",
      limits   = c(1, fill_limit),
      oob      = scales::squish,
      na.value = NA
    ) +
    scale_x_continuous(limits = ax_lim, breaks = ax_breaks) +
    scale_y_continuous(limits = ax_lim, breaks = ax_breaks) +
    coord_fixed() +
    facet_wrap(~ strip, nrow = 1) +
    labs(x = NULL, y = NULL) +
    theme_storyline()
  
  if (!show_y_text) {
    p <- p + theme(axis.text.y  = element_blank(),
                   axis.ticks.y = element_blank())
  }
  p
}

# One panel for (season, row_group, model)
make_panel <- function(season, row_group, model, strip_label, first_in_row,
                       pairs_all, scale_args) {
  df_panel <- filter(pairs_all, season == !!season,
                     row_group == !!row_group, col_label == !!model)
  build_panel_plot(
    df_panel, strip_label,
    scale_args$ax_lim, scale_args$ax_breaks,
    scale_args$common_breaks, scale_args$fill_limit,
    show_y_text = first_in_row || !hide_inner_y_labels
  )
}

# ---------------------------------------------------------------
# ONE ROW BLOCK: title on top, shared y title on the left, shared
# x title underneath, panels in the middle. Title / axis-title rows
# use fixed cm sizes so they never eat into the (square) panels.
# ---------------------------------------------------------------
build_block <- function(panels, title) {
  wrap_plots(
    list(
      label_plot(title, x = 0, hjust = 0, size = 4.4, face = "bold",
                 colour = row_title_colour),
      label_plot(y_axis_title, angle = 90),
      wrap_plots(panels, nrow = 1),
      label_plot(x_axis_title)
    ),
    design  = "AA\nBC\n#D",
    widths  = unit(c(0.5, 1),       c("cm", "null")),
    heights = unit(c(0.6, 1, 0.5),  c("cm", "null", "cm"))
  )
}

# Stack blocks; heights proportional to 1 / n_panels (square panels)
stack_blocks <- function(blocks, n_panels, legend_key_height = 1.6) {
  h <- 1 / n_panels + row_height_pad
  wrap_plots(blocks, ncol = 1, heights = unit(h, "null")) +
    plot_layout(guides = "collect") &
    theme(
      legend.position    = "right",
      legend.box.spacing = unit(0.1, "cm"),
      legend.key.height  = unit(legend_key_height, "cm"),
      legend.key.width   = unit(0.35, "cm")
    )
}

# ---------------------------------------------------------------
# FIGURE A: ANNUAL
# ---------------------------------------------------------------
make_annual_plot <- function(pairs_all, scale_args) {
  blocks <- list(); n_panels <- c()
  
  for (rg in row_levels) {
    cfg <- filter(model_config, row_group == rg)
    panels <- lapply(seq_len(nrow(cfg)), function(i) {
      make_panel("Annual", rg, cfg$model[i], cfg$model[i],
                 first_in_row = (i == 1), pairs_all, scale_args)
    })
    blocks   <- c(blocks, list(build_block(panels, rg)))
    n_panels <- c(n_panels, nrow(cfg))
  }
  stack_blocks(blocks, n_panels, legend_key_height = 1.35)
}

# ---------------------------------------------------------------
# FIGURE B: SUMMER + WINTER COMBINED
#   Storyline rows: 4 columns (Summer HCLIM, Summer RACMO,
#                              Winter HCLIM, Winter RACMO)
#   ERA5: a Summer row (3 panels) and a Winter row (3 panels)
# ---------------------------------------------------------------
make_summer_winter_plot <- function(pairs_all, scale_args) {
  blocks <- list(); n_panels <- c()
  sw <- c("Summer", "Winter")
  
  # Storyline 1 / 2: season x model across 4 columns
  for (rg in c("Storyline 1", "Storyline 2")) {
    cfg   <- filter(model_config, row_group == rg)
    combo <- expand.grid(model = cfg$model, season = sw,
                         stringsAsFactors = FALSE)
    # expand.grid varies model fastest -> Summer HCLIM, Summer RACMO,
    # Winter HCLIM, Winter RACMO
    panels <- lapply(seq_len(nrow(combo)), function(i) {
      make_panel(combo$season[i], rg, combo$model[i],
                 paste0(combo$model[i], " (", combo$season[i], ")"),
                 first_in_row = (i == 1), pairs_all, scale_args)
    })
    blocks   <- c(blocks, list(build_block(panels, rg)))
    n_panels <- c(n_panels, nrow(combo))
  }
  
  # ERA5: one row per season, 3 columns each
  cfg <- filter(model_config, row_group == "ERA5 Re-analysis")
  for (s in sw) {
    panels <- lapply(seq_len(nrow(cfg)), function(i) {
      make_panel(s, "ERA5 Re-analysis", cfg$model[i], cfg$model[i],
                 first_in_row = (i == 1), pairs_all, scale_args)
    })
    blocks   <- c(blocks, list(build_block(panels, paste0("ERA5 Re-analysis – ", s))))
    n_panels <- c(n_panels, nrow(cfg))
  }
  
  stack_blocks(blocks, n_panels)
}

# ---------------------------------------------------------------
# GLOBAL SCALE (Annual + Summer + Winter, shared by both figures so
# the fill legend means the same thing everywhere; capped at 400)
# ---------------------------------------------------------------
season_pairs_cache <- setNames(lapply(seasons, build_season_pairs), seasons)
all_pairs <- bind_rows(season_pairs_cache)

raw_lim   <- range(c(all_pairs$x, all_pairs$y), na.rm = TRUE)
ax_lim    <- c(floor(raw_lim[1] / 5) * 5, ceiling(raw_lim[2] / 5) * 5)
ax_breaks <- seq(ax_lim[1], ax_lim[2], by = 10)

common_breaks <- seq(ax_lim[1], ax_lim[2], length.out = 151)  # 150 bins
fill_limit    <- min(get_shared_fill_limit(all_pairs, common_breaks), 400)

scale_args <- list(ax_lim = ax_lim, ax_breaks = ax_breaks,
                   common_breaks = common_breaks, fill_limit = fill_limit)

# ---------------------------------------------------------------
# SAVE
# ---------------------------------------------------------------
a4_width  <- 7.8    # inches
a4_dpi    <- 320

# Annual: 2 + 2 + 3 panel rows (FIGURE 5 Main MANUSCRIPT)
ggsave(
  file.path(outpath, "AntAirICE_vs_Models_annual.png"),
  make_annual_plot(all_pairs, scale_args),
  width = 6.6, height = 9.0, dpi = a4_dpi, bg = "white"   # narrower -> legend sits closer
)
message("Annual plot saved.")

# Summer + Winter combined: 4 blocks, so a bit taller (A4 max ~11.7 in)
# FIGURE 6 MAIN MANUSCRIPT
ggsave(
  file.path(outpath, "AntAirICE_vs_Models_summer_winter.png"),
  make_summer_winter_plot(all_pairs, scale_args),
  width = a4_width, height = 10.8, dpi = a4_dpi, bg = "white"
)
message("Summer/Winter combined plot saved.")


# ###############################################################
# FIGURE 2: BIAS RASTER MAPS (Storyline 1 / 2 / ERA5 rows)
# ###############################################################

# ---------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------
rast_to_df <- function(r, val_name = "value") {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- val_name
  df
}

agg_fact <- 10  # aggregation factor (coarser cells stand out more)

# ---------------------------------------------------------------
# DIVERGING COLOUR SCALE (zero always lines up with white)
# ---------------------------------------------------------------
diverging_ramp <- c(
  "#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
  "white",
  "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"
)

make_scale_values <- function(min_val, max_val) {
  white_pos <- (0 - min_val) / (max_val - min_val)
  c(
    0,
    white_pos * 0.25, white_pos * 0.50, white_pos * 0.75, white_pos * 0.95,
    white_pos,
    white_pos + (1 - white_pos) * 0.15,
    white_pos + (1 - white_pos) * 0.40,
    white_pos + (1 - white_pos) * 0.65,
    white_pos + (1 - white_pos) * 0.85,
    1
  )
}

# Slightly larger colour bar (was 3.2 x 0.4 cm)
make_cbar <- function() {
  guide_colorbar(
    barheight       = unit(3.9, "cm"),
    barwidth        = unit(0.5, "cm"),
    ticks.colour    = "black",
    frame.colour    = "black",
    frame.linewidth = 0.4
  )
}

# ---------------------------------------------------------------
# STYLE
# ---------------------------------------------------------------
theme_bias_map <- function(base_size = 11, base_family = "Helvetica") {
  theme_void(base_size = base_size) +
    theme(
      text         = element_text(family = base_family),
      plot.title   = element_text(face = "bold", size = 10.5, hjust = 0,
                                  colour = "grey20", margin = margin(b = 4)),
      legend.title = element_text(size = 9.5, colour = "grey20"),
      legend.text  = element_text(size = 8.5, colour = "grey40"),
      plot.margin  = margin(4, 4, 4, 4)
    )
}

# ---------------------------------------------------------------
# COASTLINE
# ---------------------------------------------------------------

# Antarctic continent, from SCAR Antarctic Digital Database high resolution vector polygons of the Antarctic coastline v7.12 (Gerrish et al., 2026) 
# NOTE: Downloaded the zipped file and unzipped it manually

coast <- st_read(
  here("Data/add_coastline_medium_res_polygon_v7_10.shp"),
  quiet = TRUE
)
coast_fill <- "grey95"

# ---------------------------------------------------------------
# MODEL CONFIGURATION (row label here is "ERA5 Reanalysis")
# ---------------------------------------------------------------
row_levels <- c("Storyline 1", "Storyline 2", "ERA5 Reanalysis")

model_config <- tribble(
  ~model,   ~driving,     ~row_group,
  "HCLIM",  "MPI_ESM1",   "Storyline 1",
  "RACMO",  "MPI_ESM1",   "Storyline 1",
  "HCLIM",  "CESM2",      "Storyline 2",
  "RACMO",  "CESM2",      "Storyline 2",
  "HCLIM",  "ERA5",       "ERA5 Reanalysis",
  "RACMO",  "ERA5",       "ERA5 Reanalysis",
  "MetUM",  "ERA5",       "ERA5 Reanalysis"
) %>%
  mutate(
    col_label = model,
    folder    = paste0(model, "_", driving)
  )

# ---------------------------------------------------------------
# BUILD BIAS RASTERS FOR ONE SEASON (Model - AntAirICE, 10 km)
# ---------------------------------------------------------------
build_season_bias <- function(season) {
  ref <- rast(ref_path(season))
  
  pmap(model_config, function(model, driving, row_group, col_label, folder) {
    mod <- rast(model_path(folder, season))
    check_same_grid(ref, mod, paste(model, driving, season))
    
    bias_r    <- mod - ref
    bias_10km <- aggregate(bias_r, fact = agg_fact, fun = "mean", na.rm = TRUE)
    
    list(row_group = row_group, col_label = col_label, r = bias_10km)
  })
}

# ---------------------------------------------------------------
# SHARED DIVERGING SCALE FOR ONE SEASON
# ---------------------------------------------------------------
get_season_scale <- function(bias_list) {
  all_vals <- unlist(lapply(bias_list, function(b) values(b$r, na.rm = TRUE)))
  diff_min <- floor(min(all_vals, na.rm = TRUE))
  diff_max <- ceiling(max(all_vals, na.rm = TRUE))
  list(
    min    = diff_min,
    max    = diff_max,
    breaks = pretty(c(diff_min, diff_max), n = 6),
    values = make_scale_values(diff_min, diff_max)
  )
}

# ---------------------------------------------------------------
# BUILD ONE BIAS PANEL
# ---------------------------------------------------------------
build_bias_panel <- function(r, col_label, scale_info) {
  df <- rast_to_df(r, "bias")
  
  ggplot() +
    geom_sf(data = coast, fill = coast_fill, colour = "grey35", linewidth = 0.15) +
    geom_tile(data = df, aes(x = x, y = y, fill = bias)) +
    scale_fill_gradientn(
      colours = diverging_ramp,
      values  = scale_info$values,
      limits  = c(scale_info$min, scale_info$max),
      breaks  = scale_info$breaks,
      name    = "Bias (°C)",
      oob     = squish,
      guide   = make_cbar()
    ) +
    coord_sf(expand = FALSE) +
    labs(title = col_label) +
    theme_bias_map()
}

# ---------------------------------------------------------------
# FLAT LAYOUT HELPERS
# Nested patchworks drop their plot_annotation() titles, and the
# 3-panel row ends up left-aligned. So instead: each row title is
# its own thin text strip, and everything goes into ONE flat
# patchwork on a 6-column grid (2-panel rows: 3 columns per panel;
# 3-panel row: 2 columns per panel) -> bottom row is centred and
# all rows span the same width.
# ---------------------------------------------------------------
title_strip <- function(txt) {
  ggplot() +
    annotate("text", x = 0, y = 0.5, label = txt, hjust = 0,
             fontface = "bold", size = 4.4, colour = "grey15",
             family = "Helvetica") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# row_panels: named list, one element per row, each a list of ggplots
assemble_rows <- function(row_panels) {
  pieces  <- list()
  areas   <- list()
  heights <- list()
  
  for (i in seq_along(row_panels)) {
    panels <- row_panels[[i]]
    span   <- 6 / length(panels)
    t_row  <- 2 * i - 1   # title strip row
    p_row  <- 2 * i       # panel row
    
    pieces <- c(pieces, list(title_strip(names(row_panels)[i])))
    areas  <- c(areas,  list(area(t_row, 1, t_row, 6)))
    
    for (j in seq_along(panels)) {
      pieces <- c(pieces, list(panels[[j]]))
      areas  <- c(areas,  list(area(p_row, (j - 1) * span + 1, p_row, j * span)))
    }
    
    # Title strip: fixed height. Panel row: relative height; 2-panel rows
    # are wider so (for roughly square maps) get 1.5x the 3-panel height.
    heights <- c(heights, list(unit(0.7, "cm"),
                               unit(6 / length(panels) / 2, "null")))
  }
  
  wrap_plots(pieces) +
    plot_layout(
      design  = do.call(c, areas),
      heights = do.call(grid::unit.c, heights),
      widths  = unit(rep(1, 6), "null"),
      guides  = "collect"
    )
}

# ---------------------------------------------------------------
# BUILD FULL SEASON BIAS FIGURE
# ---------------------------------------------------------------
build_bias_row_panels <- function(row_title, bias_list, scale_info) {
  row_items <- Filter(function(b) b$row_group == row_title, bias_list)
  lapply(row_items, function(b) build_bias_panel(b$r, b$col_label, scale_info))
}

make_season_bias_plot <- function(season) {
  bias_list  <- build_season_bias(season)
  scale_info <- get_season_scale(bias_list)
  
  row_panels <- setNames(
    lapply(row_levels, build_bias_row_panels, bias_list, scale_info),
    row_levels
  )
  
  assemble_rows(row_panels) &
    theme(
      legend.position    = "right",
      legend.box.spacing = unit(0.9, "cm"),   # distance from plots
      legend.title       = element_text(size = 10.5, colour = "grey20"),
      legend.text        = element_text(size = 9.5,  colour = "grey40")
    )
}

# ---------------------------------------------------------------
# RUN FOR EACH SEASON AND SAVE (A4-friendly)
# ---------------------------------------------------------------
a4_width  <- 7.8
a4_height <- 10.5
a4_dpi    <- 320

for (season in seasons) {
  message("Processing bias maps: ", season, " ...")
  
  ggsave(
    file.path(outpath, paste0("AntAirICE_vs_Models_bias_", tolower(season), ".png")),
    make_season_bias_plot(season),
    width  = a4_width,
    height = a4_height,
    dpi    = a4_dpi,
    bg     = "white"
  )
  message(season, " bias plot saved.")
}


# ###############################################################
# PART 3: SUMMARY STATISTICS
# Mean bias, RMSE, MAE and R2
# ###############################################################

calculate_metrics <- function(ref, mod, label = "") {
  check_same_grid(ref, mod, label)
  
  vals <- data.frame(
    actual    = as.vector(values(ref)),
    predicted = as.vector(values(mod))
  )
  vals <- vals[complete.cases(vals), ]
  
  lm_fit <- lm(predicted ~ actual, data = vals)
  
  data.frame(
    n_grid_cells = nrow(vals),
    mean_bias    = mean(vals$predicted - vals$actual),   # model - AntAirICE
    RMSE         = Metrics::rmse(vals$actual, vals$predicted),
    MAE          = Metrics::mae(vals$actual, vals$predicted),
    R2           = summary(lm_fit)$r.squared
  )
}

metrics_list <- list()

for (season in seasons) {
  message("Calculating metrics for ", season, " ...")
  
  ref <- rast(ref_path(season))
  
  metrics_df <- pmap_dfr(
    model_config,
    function(model, driving, row_group, col_label, folder) {
      mod     <- rast(model_path(folder, season))
      metrics <- calculate_metrics(ref, mod, paste(model, driving, season))
      
      metrics$model     <- model
      metrics$driving   <- driving
      metrics$row_group <- row_group
      metrics$season    <- season
      metrics
    }
  ) %>%
    select(season, row_group, model, driving,
           n_grid_cells, mean_bias, RMSE, MAE, R2)
  
  write.csv(
    metrics_df,
    file.path(outpath, paste0("AntAirICE_model_metrics_", tolower(season), ".csv")),
    row.names = FALSE
  )
  metrics_list[[season]] <- metrics_df
  message(season, " metrics saved.")
}

# Combined summary: Annual + Summer + Winter in one spreadsheet
write.csv(
  bind_rows(metrics_list),
  file.path(outpath, "AntAirICE_model_metrics_all_seasons.csv"),
  row.names = FALSE
)
message("Combined metrics (Annual, Summer, Winter) saved.")


# ###############################################################
# PART 3B: MEAN BIAS HEATMAP
#
# Rows:          Storyline 1 / Storyline 2 / ERA5
# Column groups: HCLIM (ALL DJF JJA) | RACMO (ALL DJF JJA) | MetUM (ALL DJF JJA)
# ###############################################################

# ---------------------------------------------------------------
# LOAD METRICS OUTPUT
# ---------------------------------------------------------------
metric_files <- setNames(
  file.path(outpath, paste0("AntAirICE_model_metrics_", tolower(seasons), ".csv")),
  seasons
)

bias_df <- bind_rows(
  lapply(names(metric_files), function(season) {
    df <- read.csv(metric_files[[season]])
    df$season <- season
    df
  })
)

# Drop the row_group column from the CSVs - it's re-added below from
# `model_structure` (avoids row_group.x / row_group.y after the join)
bias_df <- select(bias_df, -any_of("row_group"))

# ---------------------------------------------------------------
# MODEL / STORYLINE STRUCTURE
# ---------------------------------------------------------------
model_structure <- tribble(
  ~row_group,     ~model,  ~driving,    ~model_group,
  "Storyline 1",  "HCLIM", "MPI_ESM1",  "HCLIM",
  "Storyline 1",  "RACMO", "MPI_ESM1",  "RACMO",
  "Storyline 2",  "HCLIM", "CESM2",     "HCLIM",
  "Storyline 2",  "RACMO", "CESM2",     "RACMO",
  "ERA5",         "HCLIM", "ERA5",      "HCLIM",
  "ERA5",         "RACMO", "ERA5",      "RACMO",
  "ERA5",         "MetUM", "ERA5",      "MetUM"
)

bias_df <- bias_df %>%
  mutate(
    season_label = recode(season,
                          "Annual" = "ALL",
                          "Summer" = "DJF",
                          "Winter" = "JJA")
  )

# ---------------------------------------------------------------
# COMPLETE 3 x 3 x 3 GRID
# (MetUM cells stay empty for the storylines - MetUM is ERA5-only)
# ---------------------------------------------------------------
plot_grid <- expand_grid(
  row_group    = c("Storyline 1", "Storyline 2", "ERA5"),
  model_group  = c("HCLIM", "RACMO", "MetUM"),
  season_label = c("ALL", "DJF", "JJA")
) %>%
  left_join(model_structure, by = c("row_group", "model_group")) %>%
  left_join(
    bias_df %>% select(model, driving, season_label, mean_bias),
    by = c("model", "driving", "season_label")
  )

plot_grid$row_group    <- factor(plot_grid$row_group,
                                 levels = rev(c("Storyline 1", "Storyline 2", "ERA5")))
plot_grid$model_group  <- factor(plot_grid$model_group,
                                 levels = c("HCLIM", "RACMO", "MetUM"))
plot_grid$season_label <- factor(plot_grid$season_label,
                                 levels = c("ALL", "DJF", "JJA"))

# ---------------------------------------------------------------
# SYMMETRIC BIAS SCALE (limit rounded up to nearest 0.5)
# ---------------------------------------------------------------
max_bias   <- max(abs(plot_grid$mean_bias), na.rm = TRUE)
bias_limit <- ceiling(max_bias * 2) / 2

# ---------------------------------------------------------------
# BUILD HEATMAP
# ---------------------------------------------------------------
bias_plot <- ggplot(plot_grid, aes(x = season_label, y = row_group)) +
  
  geom_tile(aes(fill = mean_bias), colour = "white", linewidth = 0.8) +
  
  geom_text(
    aes(
      label  = ifelse(is.na(mean_bias), "", sprintf("%.2f", mean_bias)),
      colour = ifelse(!is.na(mean_bias) & abs(mean_bias) > bias_limit * 0.45,
                      "white", "black")
    ),
    size = 4
  ) +
  
  scale_colour_identity() +
  
  scale_fill_gradient2(
    low      = "#2C7BB6",
    mid      = "white",
    high     = "#D7191C",
    midpoint = 0,
    limits   = c(-bias_limit, bias_limit),
    name     = "Mean bias (°C)",
    na.value = "grey95"
  ) +
  
  facet_grid(. ~ model_group) +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text.y       = element_text(size = 11.5, face = "bold", colour = "black"),
    axis.text.x       = element_text(size = 10.5, colour = "black"),
    axis.ticks        = element_line(colour = "grey60", linewidth = 0.35),
    axis.line         = element_blank(),
    strip.background  = element_blank(),
    strip.text        = element_text(size = 10.5, colour = "black"),
    panel.spacing.x   = unit(0.8, "lines"),   # gap between column groups
    legend.title      = element_text(size = 10.5),
    legend.text       = element_text(size = 9.5),
    legend.key.height = unit(1, "cm"),
    legend.key.width  = unit(0.35, "cm"),
    plot.margin       = margin(5, 5, 5, 5)
  )

# ---------------------------------------------------------------
# SAVE
# ---------------------------------------------------------------
ggsave(
  file.path(outpath, "AntAirICE_model_bias_heatmap_grouped.png"),
  bias_plot,
  width  = 7.5,
  height = 3.8,
  dpi    = 300
)

###############################################
# References --------------------------------------------------------------
################################################

# Gerrish, L., Ireland, L., Fretwell, P., Cooper, P., & Skachkova, A. (2026). High resolution vector polygons of the Antarctic coastline 2 (Version 7.12) [Data set]. NERC EDS UK Polar Data Centre. https://doi.org/10.5285/13c4d2f1-8903-4d7f-8977-592121975554 [Accessed 09 June 2026]