
# ==============================================================================
# FIGURE 4 main manuscript — Temperature difference (future − historical)
# Annual / Summer / Winter
# ==============================================================================#
# INPUTS:
# 1. Regional climate model (RCM) rasters -- ice-free land, regridded
#    Difference between future and historical rasters
#    Created with the 3_Storyline_Summaries.R script
# 2. Scope boundary (splits ice-free land into Peninsula vs Continent)
#    Polygon created by the authors.
#    Available in Data/Peninsula_Continent_Boundary.shp
#    Scopes plotted: All (ice-free land), Peninsula, Continent


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)  
library(here)
library(svglite)


# ------------------------------------------------------------------
# 0. SETTINGS
# ------------------------------------------------------------------

# Final variable set + display order (left -> right, top -> bottom)
subset_vars <- c(
  "AnnualTemp", "SummerTemp", "WinterTemp", "DegreeDays0",
  "TotalAnnualPrecip", "WindSpeed", "SeaIceConcBuffer", "RelHumidity"
)

# Nice facet labels (edit to match your variable naming / units)
var_labels <- c(
  AnnualTemp        = "Annual temp (\u00B0C)",
  SummerTemp        = "Summer temp (\u00B0C)",
  WinterTemp        = "Winter temp (\u00B0C)",
  DegreeDays0       = "Degree days <0\u00B0C",
  TotalAnnualPrecip = "Total annual precip (mm)",
  WindSpeed         = "Wind speed (m/s)",
  SeaIceConcBuffer        = "Mean sea ice conc. (%)",
  RelHumidity       = "Relative humidity (%)"
)

driver_colours <- c(CESM2 = "#eb6834", MPI_ESM1 = "#2a78d6")

# Paler tints of the driver colours, used for the individual-RCM ticks so
# each tick reads as "belonging to" its driver without a separate legend
# entry (matches the softer orange/blue used in the annotated mock-up)
pale_driver_colours <- c(CESM2 = "#f3b79a", MPI_ESM1 = "#a9c9ee")

# Neutral grey for the driver-mean-to-driver-mean bar and the full-width
# baseline. Both share one linewidth (axis_linewidth below) so the
# "shaded" (data range) and "unshaded" (rest of the axis) portions of the
# horizontal line read as one continuous track, just different colours.
rcm_bar_colour   <- "grey60"
baseline_colour  <- "grey88"
axis_linewidth   <- 3
tick_linewidth   <- 1.1

# Point + text styling
point_size    <- 4.5
point_stroke  <- 1.1
title_colour  <- "grey45"   # panel/scope titles - light grey, slightly
axis_num_colour <- "grey60" # darker than the axis numbers
label_family <- "sans"      # swap for a specific installed font if you want
# an exact match, e.g. "Helvetica Neue" / "Roboto"

# --- zero-reference line (for wind speed) --------------------
# A small, thin black vertical line marking x = 0, drawn only in the
# panel(s) named here. Thinner than the RCM ticks (tick_linewidth) so it
# reads as a subtle reference mark rather than another data element.
zero_line_vars     <- c("WindSpeed")
zero_line_colour   <- "black"
zero_line_linewidth <- 0.5   # thinner than tick_linewidth (1.1)

# --- mean-diff value labels under each driver-mean point ------
# Slightly darker than the axis number colour (axis_num_colour, grey60)
# so the exact values read as a touch more prominent, but still small -
# sized to roughly match the axis text (rel(0.85) of base_size 11).
label_colour     <- "grey35"
label_text_size  <- 2.6   # slightly smaller than the axis text (was 3.3)
label_y_pos      <- 0.78  # below the y = 1 track line

# Variables whose x-axis should be forced to a shared range so the
# panels can be visually compared to one another
shared_axis_group <- c("AnnualTemp", "SummerTemp", "WinterTemp")

# Axis-minimum "floors" -- applied after the automatic 0-floor logic
# below. The axis will extend further than this if the data needs it
# (e.g. an RCM below -0.1), but will never be *less* extreme than this
# even if all the data happens to sit above it.
axis_min_floor <- c(WindSpeed = -0.1)

# --- Axis-maximum "ceilings" -- mirror of axis_min_floor above,
# but for the upper (right-hand) limit. Forces the axis max UP to at
# least the given value even if the data itself doesn't reach it (e.g.
# SeaIceConcBuffer's diffs are all negative, with a data max of -1, but
# we still want the right edge of that panel's axis anchored at 0).
axis_max_ceiling <- c(SeaIceConcBuffer = 0)

# ------------------------------------------------------------------
# 1. DATA PREP FUNCTION (per scope)
# ------------------------------------------------------------------

prepare_storyline_data <- function(model_df, subset_vars) {
  
  present_vars <- intersect(subset_vars, unique(model_df$variable))
  missing_vars <- setdiff(subset_vars, present_vars)
  
  model_means <- model_df %>%
    filter(variable %in% present_vars) %>%
    pivot_longer(
      cols = c(HCLIM_CESM2, HCLIM_MPI_ESM1, RACMO_CESM2, RACMO_MPI_ESM1),
      names_to = "Model",
      values_to = "Value"
    ) %>%
    mutate(
      Value  = as.numeric(sub("\\s*\\+/-.*$", "", Value)),
      RCM    = case_when(grepl("^HCLIM", Model) ~ "HCLIM",
                         grepl("^RACMO", Model) ~ "RACMO"),
      Driver = case_when(grepl("CESM2$", Model)    ~ "CESM2",
                         grepl("MPI_ESM1$", Model) ~ "MPI_ESM1")
    ) %>%
    select(Model, Driver, RCM, variable, Value) %>%
    rename(Variable = variable)
  
  # Add placeholder rows for any requested variable not yet in the data
  # (e.g. SeaIceConcBuffer) so the facet panel still reserves its slot.
  if (length(missing_vars) > 0) {
    placeholder <- expand.grid(
      Model = NA_character_, Driver = c("CESM2", "MPI_ESM1"),
      RCM = NA_character_, Variable = missing_vars,
      stringsAsFactors = FALSE
    ) %>% mutate(Value = NA_real_)
    model_means <- bind_rows(model_means, placeholder)
  }
  
  model_means <- model_means %>%
    mutate(Variable = factor(Variable, levels = subset_vars))
  
  storyline_means <- model_means %>%
    group_by(Driver, Variable) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(Variable = factor(Variable, levels = subset_vars))
  
  segment_df <- storyline_means %>%
    pivot_wider(names_from = Driver, values_from = Value) %>%
    mutate(Variable = factor(Variable, levels = subset_vars))
  
  list(
    model_means     = model_means,
    storyline_means = storyline_means,
    segment_df      = segment_df,
    missing_vars    = missing_vars
  )
}

# ------------------------------------------------------------------
# 2. PER-VARIABLE X-AXIS SCALES (0 fixed at left, single max break,
#    shared range for the temperature trio, manual overrides supported)
# ------------------------------------------------------------------

build_x_scales <- function(model_means, subset_vars, shared_axis_group,
                           axis_min_floor = c(), axis_max_ceiling = c()) {
  
  axis_range <- model_means %>%
    group_by(Variable) %>%
    summarise(min_val = min(Value, na.rm = TRUE),
              max_val = max(Value, na.rm = TRUE))
  
  axis_min <- tibble::deframe(axis_range %>% select(Variable, min_val))
  axis_max <- tibble::deframe(axis_range %>% select(Variable, max_val))
  
  # replace any +/-Inf (all-NA placeholder variables) with nominal values
  axis_min[!is.finite(axis_min)] <- 0
  axis_max[!is.finite(axis_max)] <- 1
  
  # 0 is always kept as the *lower bound* unless the data actually goes
  # negative (e.g. WindSpeed), in which case drop the lower limit to fit it
  axis_min <- pmin(axis_min, 0)
  
  # apply any floors (e.g. WindSpeed -> at least -0.1) -- takes whichever
  # is more extreme so a lower RCM value still gets fully captured
  for (v in names(axis_min_floor)) {
    if (v %in% names(axis_min)) {
      axis_min[[v]] <- min(axis_min[[v]], axis_min_floor[[v]])
    }
  }
  
  # apply any ceilings (e.g. SeaIceConcBuffer -> right edge anchored at
  # 0 even though the data itself is all negative) -- takes whichever is
  # LESS extreme (i.e. pulls the max UP towards/ to the ceiling value)
  for (v in names(axis_max_ceiling)) {
    if (v %in% names(axis_max)) {
      axis_max[[v]] <- max(axis_max[[v]], axis_max_ceiling[[v]])
    }
  }
  
  shared_max <- max(axis_max[shared_axis_group], na.rm = TRUE)
  axis_max[shared_axis_group] <- shared_max
  
  range_df <- tibble::tibble(
    Variable = factor(subset_vars, levels = subset_vars),
    xmin     = axis_min[subset_vars],
    xmax     = axis_max[subset_vars]
  )
  
  # Build the scale list IN THE SAME ORDER as subset_vars/facet levels
  scales <- lapply(subset_vars, function(v) {
    lower <- axis_min[[v]]
    upper <- axis_max[[v]]
    scale_x_continuous(
      limits = c(lower, upper),
      breaks = c(round(lower, 1), round(upper, 1)),
      expand = expansion(mult = c(0.03, 0.08))
    )
  })
  
  list(scales = scales, range_df = range_df)
}

# ------------------------------------------------------------------
# 3. PLOTTING FUNCTION
# ------------------------------------------------------------------

make_storyline_plot <- function(model_df, subset_vars, var_labels,
                                driver_colours, shared_axis_group,
                                axis_min_floor = c(), axis_max_ceiling = c(),
                                zero_line_vars = c(), scope_title = NULL) {
  
  d  <- prepare_storyline_data(model_df, subset_vars)
  xs <- build_x_scales(d$model_means, subset_vars, shared_axis_group,
                       axis_min_floor, axis_max_ceiling)
  
  # One row per panel that should get a zero-reference line
  zero_line_df <- tibble::tibble(
    Variable = factor(intersect(zero_line_vars, subset_vars), levels = subset_vars),
    x        = 0
  )
  
  p <- ggplot() +
    
    # pale full-width baseline ("unshaded" track), one per panel, spanning
    # that panel's axis -- same linewidth as the shaded bar below so the
    # two read as one continuous line
    geom_segment(data = xs$range_df,
                 aes(x = xmin, xend = xmax, y = 1, yend = 1),
                 colour = baseline_colour, linewidth = axis_linewidth,
                 lineend = "round") +
    
    # grey bar ("shaded" portion) joining the two driver (storyline) means
    geom_segment(data = d$segment_df,
                 aes(x = CESM2, xend = MPI_ESM1, y = 1, yend = 1),
                 colour = rcm_bar_colour, linewidth = axis_linewidth,
                 lineend = "round") +
    
    # NEW: thin black zero-reference line (e.g. wind speed), drawn under
    # the RCM ticks/points so it reads as a subtle background reference
    geom_segment(data = zero_line_df,
                 aes(x = x, xend = x, y = 0.8, yend = 1.2),
                 colour = zero_line_colour, linewidth = zero_line_linewidth,
                 na.rm = TRUE) +
    
    # individual RCM marks as short vertical ticks, coloured by driver
    # (pale blue for CESM2-driven runs, pale orange for MPI_ESM1-driven)
    geom_segment(data = d$model_means,
                 aes(x = Value, xend = Value, y = 0.9, yend = 1.1,
                     colour = Driver),
                 linewidth = tick_linewidth, na.rm = TRUE, show.legend = FALSE) +
    
    # driver (storyline) means as larger dots with a white outline
    geom_point(data = d$storyline_means,
               aes(x = Value, y = 1, fill = Driver),
               shape = 21, colour = "white", stroke = point_stroke,
               size = point_size, na.rm = TRUE) +
    
    # NEW: exact mean-diff value printed under each driver-mean point
    geom_text(data = d$storyline_means,
              aes(x = Value, label = sprintf("%.1f", Value)),
              y = label_y_pos, colour = label_colour, size = label_text_size,
              family = label_family, na.rm = TRUE) +
    
    facet_wrap2(~ Variable, scales = "free_x", ncol = 4, drop = FALSE,
                labeller = as_labeller(var_labels)) +
    facetted_pos_scales(x = xs$scales) +
    
    scale_colour_manual(values = pale_driver_colours, guide = "none") +
    scale_fill_manual(values = driver_colours, name = NULL) +
    guides(fill = guide_legend(override.aes = list(size = point_size))) +
    
    scale_y_continuous(limits = c(0.5, 1.5), breaks = NULL) +
    labs(x = NULL, y = NULL, title = scope_title) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid           = element_blank(),
      strip.text           = element_text(face = "plain", hjust = 0,
                                          colour = title_colour,
                                          family = label_family,
                                          size = rel(1.05),
                                          margin = margin(b = 3)),
      strip.placement       = "outside",
      axis.text.x           = element_text(colour = axis_num_colour,
                                           family = label_family,
                                           size = rel(0.85)),
      axis.ticks.x          = element_blank(),
      axis.line.x           = element_blank(),
      panel.spacing.x       = unit(2.2, "lines"),
      panel.spacing.y       = unit(2, "lines"),
      legend.position        = "top",
      legend.justification   = "left",
      legend.text             = element_text(family = label_family, size = rel(0.9)),
      plot.title               = element_text(face = "plain", colour = title_colour,
                                              family = label_family),
      plot.title.position      = "plot"
    )
  
  # Placeholder label for any variable not yet in the data (e.g. sea ice)
  if (length(d$missing_vars) > 0) {
    placeholder_range <- xs$range_df %>% filter(Variable %in% d$missing_vars)
    placeholder_df <- placeholder_range %>%
      mutate(x = (xmin + xmax) / 2, y = 1)
    p <- p + geom_text(data = placeholder_df,
                       aes(x = x, y = y, label = "Data pending"),
                       colour = "grey60", size = 3, fontface = "italic")
  }
  
  p
}

# ------------------------------------------------------------------
# 4. RUN FOR THE THREE SCOPES
# ------------------------------------------------------------------

all_model_df       <- read.csv(here("Data/Environmental_predictors/PolarRes26/Scenario_summaries/FUTURE_DIFF_All_by_model.csv"))
peninsula_model_df  <- read.csv(here("Data/Environmental_predictors/PolarRes26/Scenario_summaries/FUTURE_DIFF_PENINSULA_by_model.csv"))
continent_model_df  <- read.csv(here("Data/Environmental_predictors/PolarRes26/Scenario_summaries/FUTURE_DIFF_CONTINENT_by_model.csv"))

# Each scope gets its own axis ranges (computed from its own data only).
# Within a given scope's plot, all 8 panels are still the same physical
# width regardless of their individual data ranges -- that's inherent to
# facet_wrap's grid layout, not something the axis ranges affect.
plot_all <- make_storyline_plot(all_model_df, subset_vars, var_labels,
                                driver_colours, shared_axis_group,
                                axis_min_floor, axis_max_ceiling,
                                zero_line_vars, scope_title = "Ice-free land")

plot_peninsula <- make_storyline_plot(peninsula_model_df, subset_vars, var_labels,
                                      driver_colours, shared_axis_group,
                                      axis_min_floor, axis_max_ceiling,
                                      zero_line_vars, scope_title = "Peninsula")

plot_continent <- make_storyline_plot(continent_model_df, subset_vars, var_labels,
                                      driver_colours, shared_axis_group,
                                      axis_min_floor, axis_max_ceiling,
                                      zero_line_vars, scope_title = "Continent")

# ------------------------------------------------------------------
# 5. EXPORT AS TRUE VECTOR SVG (identical size => panels stay aligned
#    when the three files are combined in Inkscape)
# ------------------------------------------------------------------

panel_width  <- 10   # inches
panel_height <- 5    # inches

ggsave(here("Plots/storyline_all.svg"),       plot_all,       device = "svg",
       width = panel_width, height = panel_height)
ggsave(here("Plots/storyline_peninsula.svg"), plot_peninsula, device = "svg",
       width = panel_width, height = panel_height)
ggsave(here("Plots/storyline_continent.svg"), plot_continent, device = "svg",
       width = panel_width, height = panel_height)