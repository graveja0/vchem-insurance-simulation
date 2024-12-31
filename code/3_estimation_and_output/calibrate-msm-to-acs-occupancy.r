# This file will take the output from either smoothed_rates_spline or 
# smoothed_rates_gp, as well as the calibration targets located in
# here::here("_data/processed/acs-calibration-targets/acs-calibration-targets_2012.rds") (pre)
# and here::here("_data/processed/acs-calibration-targets/acs-calibration-targets_2022.rds") (post)
# to calibrate parameters/hyper-parameters using Incremental Mixture Importance sampling (IMIS)
if (!exists("manifest")) 
    source(here::here("code/0_setup/manifest.r"))
if (!exists("params"))
    source(here::here("code/0_setup/define-parameters.r"))

transition_rates <- read_rds(here::here("results/model_objects/transition_rates.rds"))

# The first step is to define which of the smoothing funcitons to use

smoothing = "spline"

# Define the smoothing functions 

# Gaussian Process smoothing
#
# Strengths:
# - Handles non-linear relationships well without imposing a specific functional form
# - Provides uncertainty estimates naturally through probabilistic framework
# - Good at interpolation between observed points
# - Robust to noise in the data
# - Can capture complex patterns in transition rates across ages
# - Can be calibrated by adjusting hyperparameters to match known totals
#
# Limitations:
# - Computationally intensive, especially with larger datasets
# - Sensitive to hyperparameter choices (kernel, sigma, variance)
# - May oversmooth sharp transitions that actually exist in the data
# - Extrapolation beyond observed age range can be unreliable
# - Assumes smoothness in transition rates which may not always hold
# - Memory intensive for large datasets due to kernel matrix calculations
# - Calibration requires careful tuning and may be difficult to achieve exact matches
# - No built-in constraints to ensure predictions match marginal totals
# - May need post-processing adjustments to align with external data
smooth_rates_gp <- function(transition_rates) {
    # Split data into groups
    grouped_data <- transition_rates %>%
        group_by(from, to, model) %>%
        group_split()
    
    # Function to smooth one group
    smooth_group <- function(group_df) {
        # Create prediction grid
        age_grid <- seq(min(group_df$age), max(group_df$age), by = 0.5)
        
        # Fit GP and predict
        gp <- gausspr(
            rate ~ age,
            data = group_df,
            kernel = "rbfdot",
            kpar = list(sigma = 2),
            var = 0.1
        )
        
        # Return smoothed predictions
        tibble(
            age = age_grid,
            rate = predict(gp, data.frame(age = age_grid)),
            from = first(group_df$from),
            to = first(group_df$to),
            model = first(group_df$model),
            transition_type = first(group_df$transition_type)
        )
    }
    
    # Map smooth_group over all groups and combine results
    smoothed_rates <- grouped_data %>%
        map_dfr(possibly(smooth_group, otherwise = NULL))
    
    return(smoothed_rates)
}

# Spline smoothing
#
# Strengths:
# - Computationally efficient compared to GP methods
# - Well-established statistical properties and theory
# - Single tuning parameter (k) makes calibration straightforward
# - Good balance of flexibility and stability
# - Built-in cross-validation for optimal smoothing
# - Works well with both sparse and dense data
# - Memory efficient even with large datasets
#
# Limitations:  
# - May struggle with very sharp transitions
# - Less flexible than GP in capturing complex patterns
# - Choice of k can still require some trial and error
# - No probabilistic interpretation unlike GP
# - Boundary effects can be problematic
# - Extrapolation beyond data range is unreliable
# - No direct way to incorporate uncertainty estimates
# - May need post-processing to ensure non-negative rates
#
# Calibration approach:
# - Start with default k=10 and adjust based on visual inspection
# - Use GCV score from gam() to guide k selection
# - Compare against known marginal rates if available
# - Cross-validate on held-out transitions
# - Consider different spline bases (e.g. P-splines) if needed
smooth_rates_spline <- function(transition_rates, k = 10) {
    # Split data into groups
    grouped_data <- transition_rates %>%
        group_by(from, to, model) %>%
        group_split()
    
    # Function to smooth one group
    smooth_group <- function(group_df) {
        # Create prediction grid
        age_grid <- seq(min(group_df$age), max(group_df$age), by = 0.5)
        
        # Fit GAM with cubic splines
        gam_fit <- gam(
            rate ~ s(age, bs = "cs", k = k),  # k controls smoothness
            data = group_df,
            method = "REML"  # Robust smoothing parameter selection
        )
        
        # Return smoothed predictions
        tibble(
            age = age_grid,
            rate = predict(gam_fit, data.frame(age = age_grid)),
            from = first(group_df$from),
            to = first(group_df$to),
            model = first(group_df$model),
            transition_type = first(group_df$transition_type)
        )
    }
    
    # Map smooth_group over all groups and combine results
    smoothed_rates <- grouped_data %>%
        map_dfr(possibly(smooth_group, otherwise = NULL))
    
    return(smoothed_rates)
}

# Apply both smoothing methods
smoothed_rates_gp <- smooth_rates_gp(transition_rates)
smoothed_rates_spline <- smooth_rates_spline(transition_rates)

# Plot comparison of all methods
plot_comparison <- function(rates_orig, rates_gp, rates_spline) {
    # Combine datasets
    rates_orig$method <- "Original"
    rates_gp$method <- "Gaussian Process"
    rates_spline$method <- "Spline"
    
    combined_rates <- bind_rows(
        rates_orig,
        rates_gp,
        rates_spline
    )
    
    ggplot(combined_rates, 
           aes(x = age, y = rate, color = to)) +
        geom_line(data = subset(combined_rates, method == "Original"), alpha = 0.5) +
        geom_line(data = subset(combined_rates, method == "Gaussian Process"), 
                  size = 1, alpha = 0.6) +
        geom_line(data = subset(combined_rates, method == "Spline"), 
                  size = 1, alpha = 0.6) +
        facet_wrap(model~from, scales = "free_y") +
        theme_minimal() +
        labs(
            title = "Transition Rates by Age and Initial Coverage Type",
            subtitle = "Comparison of Original, GP, and Spline Smoothing",
            x = "Age",
            y = "Transition Rate",
            color = "To Coverage",
            linetype = "Period"
        ) +
        theme(
            legend.position = "bottom",
            strip.background = element_rect(fill = "lightgray"),
            strip.text = element_text(face = "bold")
        )
}

# Create transtiion intensity plot
plot_comparison(transition_rates, smoothed_rates_gp, smoothed_rates_spline)





if (smoothing == "spline") rates <- smoothed_rates_spline 
if (smoothing == "gaussian") rates <- smoothed_rates_gp

