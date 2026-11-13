# ==============================================================================
# FIGURE 9 main manuscript — Temperature difference (future − historical)
# Annual / Summer / Winter, with CHELSA added as a 5-GCM reference
# ==============================================================================
# Differs from the FIGURE 4 main manuscript script in three ways:
#   1. Only the three temperature-difference variables are plotted.
#   2. HCLIM/RACMO values are computed directly from the "comparison" period
#      rasters (FUTURE 2071_2100 vs HISTORICAL 1985_2010), NOT from the
#      pre-tabulated FUTURE_DIFF csvs (which used 2081_2100) -- this keeps
#      the RCM period aligned with the CHELSA comparison period.
#   3. CHELSA is added as its own "storyline": each of its 5 GCM-driven
#      future rasters is diffed against the single CHELSA historical
#      raster, masked to the same scope, and plotted as pale purple ticks;
#      the mean of those 5 diff rasters is plotted as a solid purple dot,
#      alongside (not connected to) the CESM2/MPI_ESM1 dumbbell.
#
# INPUTS:
# 1. Regional climate model (RCM) rasters -- ice-free land, regridded (from script 2c)
#    Data/Environmental_predictors/PolarRes26/Regridded/<folder>/comparison/
#      <stem>_FUTURE_2071_2100_ICEFREE.tif       (future,     2071-2100)
#      <stem>_HISTORICAL_1985_2010_ICEFREE.tif   (historical, 1985-2010)
#    <folder> = HCLIM_CESM2, HCLIM_MPI_ESM1, RACMO_CESM2, RACMO_MPI_ESM1
#    (2 RCMs x 2 driving GCMs = 4 folders; 3 variables x 2 periods each)

# 2. CHELSA rasters -- ice-free land (from script 2c)
#    Data/CHELSA/comparison/
#      <stem>_FUTURE_<gcm>_2071_2100_ICEFREE.tif (future,     2071-2100)
#      <stem>_HISTORICAL_1981_2010_ICEFREE.tif   (historical, 1981-2010;
#                                                 single raster shared by all GCMs)
#    <gcm> = gfdl-esm4, ipsl-cm6a-lr, mpi-esm1-2-hr, mri-esm2-0, ukesm1-0-ll
#
# 3. Scope boundary (splits ice-free land into Peninsula vs Continent)
#    Polygon created by the authors.
#    Data/Peninsula_Continent_Boundary.shp
#    Scopes plotted: All (ice-free land), Peninsula, Continent


library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(ggh4x)
library(terra)
library(here)
library(svglite)

# ------------------------------------------------------------------
# 0. SETTINGS
# ------------------------------------------------------------------

subset_vars <- c("AnnualTemp", "SummerTemp", "WinterTemp")

var_labels <- c(
  AnnualTemp = "Annual temp diff (\u00B0C)",
  SummerTemp = "Summer temp diff (\u00B0C)",
  WinterTemp = "Winter temp diff (\u00B0C)"
)

# Filename stems (season -> raster variable name), shared by RCM + CHELSA
var_file_stem <- c(
  AnnualTemp = "Mean_Annual_Temperature",
  SummerTemp = "Mean_Summer_Temperature",
  WinterTemp = "Mean_Winter_Temperature"
)

rcm_folders <- c("HCLIM_CESM2", "HCLIM_MPI_ESM1", "RACMO_CESM2", "RACMO_MPI_ESM1")
rcm_future_range   <- "2071_2100"
rcm_hist_range     <- "1985_2010"

chelsa_gcms       <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr",
                       "mri-esm2-0", "ukesm1-0-ll")
chelsa_future_range <- "2071_2100"
chelsa_hist_range   <- "1981_2010"

# All three temp variables share one axis range/breaks so the panels are
# visually comparable (as in the original script)
shared_axis_group <- subset_vars

driver_colours      <- c(MPI_ESM1 = "#2a78d6", CESM2 = "#eb6834", CHELSA = "#7b3294")
pale_driver_colours <- c(MPI_ESM1 = "#a9c9ee", CESM2 = "#f3b79a", CHELSA = "#d0b3e0")

rcm_bar_colour   <- "grey60"
baseline_colour  <- "grey88"
axis_linewidth   <- 3
tick_linewidth   <- 1.1

point_size      <- 4.5
point_stroke    <- 1.1
title_colour    <- "grey45"
axis_num_colour <- "grey60"
label_family    <- "sans"

boundary_path <- here("Data/Peninsula_Continent_Boundary.shp")
scopes <- c("All", "Peninsula", "Continent")

# ------------------------------------------------------------------
# 1. PATH BUILDERS
# ------------------------------------------------------------------

rcm_path <- function(folder, variable, period, range_label) {
  here("Data/Environmental_predictors/PolarRes26/Regridded", folder, "comparison",
       paste0(var_file_stem[[variable]], "_", period, "_", range_label, "_ICEFREE.tif"))
}

chelsa_future_path <- function(variable, gcm) {
  here("Data/CHELSA/comparison",
       paste0(var_file_stem[[variable]], "_FUTURE_", gcm, "_", chelsa_future_range, "_ICEFREE.tif"))
}

chelsa_hist_path <- function(variable) {
  here("Data/CHELSA/comparison",
       paste0(var_file_stem[[variable]], "_HISTORICAL_", chelsa_hist_range, "_ICEFREE.tif"))
}

read_raster_safe <- function(path) {
  if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
  rast(path)
}

# Resample b onto a's grid if they don't already line up
align_rasters <- function(a, b) {
  if (!compareGeom(a, b, stopOnError = FALSE)) b <- resample(b, a, method = "bilinear")
  b
}

# ------------------------------------------------------------------
# 2. BUILD FULL-DOMAIN DIFF RASTERS (future - historical), UNMASKED
#    -- done once, then masked per-scope below, so rasters are only
#    read from disk a single time each
# ------------------------------------------------------------------

message("Loading RCM rasters and computing diffs...")
rcm_diff_rasters <- list()   # [[folder]][[variable]] -> SpatRaster

for (folder in rcm_folders) {
  for (v in subset_vars) {
    fut  <- read_raster_safe(rcm_path(folder, v, "FUTURE", rcm_future_range))
    hist <- read_raster_safe(rcm_path(folder, v, "HISTORICAL", rcm_hist_range))
    if (is.null(fut) || is.null(hist)) next
    hist <- align_rasters(fut, hist)
    rcm_diff_rasters[[folder]][[v]] <- fut - hist
  }
}

message("Loading CHELSA rasters and computing diffs...")
chelsa_diff_rasters <- list()        # [[gcm]][[variable]] -> SpatRaster
chelsa_ensemble_diff <- list()       # [[variable]]        -> SpatRaster (5-GCM mean)

for (v in subset_vars) {
  hist <- read_raster_safe(chelsa_hist_path(v))
  if (is.null(hist)) next
  
  gcm_diffs <- list()
  for (gcm in chelsa_gcms) {
    fut <- read_raster_safe(chelsa_future_path(v, gcm))
    if (is.null(fut)) next
    hist_aligned <- align_rasters(fut, hist)
    d <- fut - hist_aligned
    chelsa_diff_rasters[[gcm]][[v]] <- d
    gcm_diffs[[gcm]] <- d
  }
  
  if (length(gcm_diffs) > 0) {
    chelsa_ensemble_diff[[v]] <- app(rast(unname(gcm_diffs)), mean, na.rm = TRUE)
  }
}

# ------------------------------------------------------------------
# 3. SCOPE MASKING
# ------------------------------------------------------------------

boundary_vect <- vect(boundary_path)

apply_scope_mask <- function(r, scope) {
  if (scope == "All") return(r)
  if (scope == "Peninsula") {
    r <- mask(r, boundary_vect)
    r <- crop(r, ext(boundary_vect))
    return(r)
  }
  if (scope == "Continent") {
    return(mask(r, boundary_vect, inverse = TRUE))
  }
  stop("Unknown scope: ", scope)
}

raster_mean_scalar <- function(r) {
  as.numeric(global(r, "mean", na.rm = TRUE)[1, 1])
}

# ------------------------------------------------------------------
# 4. ASSEMBLE PLOT DATA FOR ONE SCOPE
# ------------------------------------------------------------------

prepare_scope_data <- function(scope) {
  
  # -- RCM ticks (one row per folder x variable) -----------------------
  model_means <- map_dfr(rcm_folders, function(folder) {
    map_dfr(subset_vars, function(v) {
      r <- rcm_diff_rasters[[folder]][[v]]
      if (is.null(r)) return(NULL)
      val <- raster_mean_scalar(apply_scope_mask(r, scope))
      tibble(
        Model    = folder,
        RCM      = str_split_i(folder, "_", 1),
        Driver   = str_remove(folder, paste0("^", str_split_i(folder, "_", 1), "_")),
        Variable = v,
        Value    = val
      )
    })
  }) %>%
    mutate(Variable = factor(Variable, levels = subset_vars),
           Driver   = factor(Driver, levels = c("CESM2", "MPI_ESM1", "CHELSA")))
  
  # -- RCM driver (storyline) means, and the dumbbell connecting them --
  storyline_means <- model_means %>%
    group_by(Driver, Variable) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  segment_df <- storyline_means %>%
    filter(Driver %in% c("CESM2", "MPI_ESM1")) %>%
    pivot_wider(names_from = Driver, values_from = Value) %>%
    mutate(Variable = factor(Variable, levels = subset_vars))
  
  # -- CHELSA ticks (one row per GCM x variable) ------------------------
  chelsa_means <- map_dfr(names(chelsa_diff_rasters), function(gcm) {
    map_dfr(subset_vars, function(v) {
      r <- chelsa_diff_rasters[[gcm]][[v]]
      if (is.null(r)) return(NULL)
      val <- raster_mean_scalar(apply_scope_mask(r, scope))
      tibble(Model = paste0("CHELSA_", gcm), RCM = "CHELSA", Driver = "CHELSA",
             Variable = v, Value = val)
    })
  }) %>%
    mutate(Variable = factor(Variable, levels = subset_vars),
           Driver   = factor(Driver, levels = c("CESM2", "MPI_ESM1", "CHELSA")))
  
  # -- CHELSA multi-model mean (purple circle) --------------------------
  chelsa_mean_df <- map_dfr(subset_vars, function(v) {
    r <- chelsa_ensemble_diff[[v]]
    if (is.null(r)) return(NULL)
    val <- raster_mean_scalar(apply_scope_mask(r, scope))
    tibble(Driver = "CHELSA", Variable = v, Value = val)
  }) %>%
    mutate(Variable = factor(Variable, levels = subset_vars),
           Driver   = factor(Driver, levels = c("CESM2", "MPI_ESM1", "CHELSA")))
  
  # combined dot layer (RCM driver means + CHELSA mean)
  dot_df <- bind_rows(storyline_means, chelsa_mean_df)
  
  # combined tick layer (RCM ticks + CHELSA ticks)
  tick_df <- bind_rows(model_means, chelsa_means)
  
  # combined values, for computing the shared axis range
  all_values_df <- bind_rows(tick_df, dot_df) %>% select(Variable, Value)
  
  list(
    tick_df       = tick_df,
    dot_df        = dot_df,
    segment_df    = segment_df,
    all_values_df = all_values_df
  )
}

# ------------------------------------------------------------------
# 5. PER-VARIABLE X-AXIS SCALES (0 fixed at left, shared range across
#    the temp trio, single max break)
# ------------------------------------------------------------------

build_x_scales <- function(all_values_df, subset_vars, shared_axis_group) {
  
  axis_range <- all_values_df %>%
    group_by(Variable) %>%
    summarise(min_val = min(Value, na.rm = TRUE),
              max_val = max(Value, na.rm = TRUE))
  
  axis_min <- tibble::deframe(axis_range %>% select(Variable, min_val))
  axis_max <- tibble::deframe(axis_range %>% select(Variable, max_val))
  
  axis_min[!is.finite(axis_min)] <- 0
  axis_max[!is.finite(axis_max)] <- 1
  
  # keep 0 in view -- extend down if any diff is negative, extend up if
  # any diff is positive (these are anomalies, so 0 is the reference)
  axis_min <- pmin(axis_min, 0)
  axis_max <- pmax(axis_max, 0)
  
  shared_min <- min(axis_min[shared_axis_group], na.rm = TRUE)
  shared_max <- max(axis_max[shared_axis_group], na.rm = TRUE)
  axis_min[shared_axis_group] <- shared_min
  axis_max[shared_axis_group] <- shared_max
  
  range_df <- tibble::tibble(
    Variable = factor(subset_vars, levels = subset_vars),
    xmin     = axis_min[subset_vars],
    xmax     = axis_max[subset_vars]
  )
  
  scales <- lapply(subset_vars, function(v) {
    lower <- axis_min[[v]]
    upper <- axis_max[[v]]
    scale_x_continuous(
      limits = c(lower, upper),
      breaks = c(round(lower, 1), 0, round(upper, 1)) %>% unique(),
      expand = expansion(mult = c(0.05, 0.08))
    )
  })
  
  list(scales = scales, range_df = range_df)
}

# ------------------------------------------------------------------
# 6. PLOTTING FUNCTION
# ------------------------------------------------------------------

make_storyline_plot <- function(scope, scope_title) {
  
  d  <- prepare_scope_data(scope)
  xs <- build_x_scales(d$all_values_df, subset_vars, shared_axis_group)
  
  ggplot() +
    
    geom_segment(data = xs$range_df,
                 aes(x = xmin, xend = xmax, y = 1, yend = 1),
                 colour = baseline_colour, linewidth = axis_linewidth,
                 lineend = "round") +
    
    geom_segment(data = d$segment_df,
                 aes(x = CESM2, xend = MPI_ESM1, y = 1, yend = 1),
                 colour = rcm_bar_colour, linewidth = axis_linewidth,
                 lineend = "round") +
    
    geom_segment(data = d$tick_df,
                 aes(x = Value, xend = Value, y = 0.95, yend = 1.05, colour = Driver),
                 linewidth = tick_linewidth, na.rm = TRUE, show.legend = FALSE) +
    
    geom_point(data = d$dot_df,
               aes(x = Value, y = 1, fill = Driver),
               shape = 21, colour = "white", stroke = point_stroke,
               size = point_size, na.rm = TRUE) +
    
    facet_wrap2(~ Variable, scales = "free_x", ncol = 3, drop = FALSE,
                labeller = as_labeller(var_labels)) +
    facetted_pos_scales(x = xs$scales) +
    
    scale_colour_manual(values = pale_driver_colours, guide = "none") +
    scale_fill_manual(values = driver_colours, name = NULL,
                      breaks = c("CESM2", "MPI_ESM1", "CHELSA")) +
    guides(fill = guide_legend(override.aes = list(size = point_size))) +
    
    scale_y_continuous(limits = c(0.5, 1.5), breaks = NULL) +
    labs(x = NULL, y = NULL, title = scope_title,
         subtitle = "Difference: future (2071\u20132100) \u2212 historical") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid           = element_blank(),
      strip.text           = element_text(face = "plain", hjust = 0,
                                          colour = title_colour,
                                          family = label_family,
                                          size = rel(1.05),
                                          margin = margin(b = 3)),
      strip.placement      = "outside",
      axis.text.x          = element_text(colour = axis_num_colour,
                                          family = label_family,
                                          size = rel(0.85)),
      axis.ticks.x         = element_blank(),
      axis.line.x          = element_blank(),
      panel.spacing.x      = unit(2.2, "lines"),
      panel.spacing.y      = unit(2, "lines"),
      legend.position      = "top",
      legend.justification = "left",
      legend.text          = element_text(family = label_family, size = rel(0.9)),
      plot.title            = element_text(face = "plain", colour = title_colour,
                                           family = label_family),
      plot.subtitle         = element_text(colour = "grey55", family = label_family,
                                           size = rel(0.85)),
      plot.title.position   = "plot"
    )
}

# ------------------------------------------------------------------
# 7. RUN FOR THE THREE SCOPES AND EXPORT
# ------------------------------------------------------------------

panel_width  <- 7
panel_height <- 4.5

scope_titles <- c(All = "Ice-free land", Peninsula = "Peninsula", Continent = "Continent")

for (scope in scopes) {
  p <- make_storyline_plot(scope, scope_titles[[scope]])
  out_file <- here("Plots", paste0("storyline_tempdiff_", tolower(scope), ".svg"))
  ggsave(out_file, p, device = "svg", width = panel_width, height = panel_height)
  message("Saved: ", out_file)
}