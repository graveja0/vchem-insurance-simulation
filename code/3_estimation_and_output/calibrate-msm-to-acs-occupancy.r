if (!exists("manifest")) 
    source(here::here("code/0_setup/manifest.r"))
if (!exists("params"))
    source(here::here("code/0_setup/define-parameters.r"))


# Calibration Targets from ACS

acs <- read_rds(
    here::here(
        "_data/processed/acs-calibration-targets/acs-calibration-targets_2012.rds"
    )
) %>%
    mutate(time = "pre") %>%
    bind_rows(read_rds(
        here::here(
            "_data/processed/acs-calibration-targets/acs-calibration-targets_2022.rds"
        )
    ) %>% mutate(time = "post")) %>% 
    mutate(type = factor(type, labels = params$v_tr_names)) %>% 
    rename(age = agep)


# Estimated Transtiion Rates

transition_rates <-  # created in examine-model-fit.r
    read_rds(here::here("results/model_objects/transition_rates.rds"))

# Plot transition rates
ggplot(transition_rates, 
       aes(x = age, y = rate, color = to, linetype = model)) +
    geom_line() +
    facet_wrap(~from, scales = "free_y") +
    theme_minimal() +
    labs(
        title = "Transition Rates by Age and Initial Coverage Type",
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


# Smooth transition rates using age-based splines

grouped_data <- 
    transition_rates %>%
    group_by(from, to, model) %>%
    group_split()

k = 20

smooth_group_spline <- function(group_df) {
    # Create prediction grid
    age_grid <- seq(min(group_df$age), max(group_df$age), by = 0.5)
    
    # Define key policy/life transition ages with more precision
    key_ages <- list(
        public_loss = c(18, 19),     # Medicaid/CHIP age-out
        college = c(22, 23),         # College graduation
        aca_young = c(25, 26),       # ACA young adult provision
        early_career = c(27, 28)     # Early career transitions
    )
    
    # Set up smoothing parameters based on transition type
    if(group_df$from[1] == "Public") {
        # Public coverage loss needs sharp transitions at eligibility boundaries
        k_value <- 30  # Higher k to allow more flexibility
        bs_type <- "ad"  # Adaptive smoothing
        
    } else if(group_df$from[1] == "Employer" || group_df$to[1] == "Employer") {
        # Employer transitions need to capture graduation and early career
        k_value <- 25
        bs_type <- "cr"  # Cubic regression spline
        
    } else if(group_df$to[1] == "Public") {
        # Transitions into public need to respect eligibility
        k_value <- 20
        bs_type <- "ad"  # Adaptive smoothing
        
    } else {
        # Default smoothing for other transitions
        k_value <- k  # Use global k value
        bs_type <- "cs"  # Cubic spline
    }
    
    # Fit GAM with appropriate basis and k value
    gam_fit <- gam(
        rate ~ s(age, bs = bs_type, k = k_value),
        data = group_df,
        method = "REML"
    )
    
    # Get initial predictions
    predictions <- predict(gam_fit, newdata = data.frame(age = age_grid))
    
    # For public and employer transitions, ensure sharp changes at key ages
    if(group_df$from[1] == "Public" || group_df$to[1] == "Public" ||
       group_df$from[1] == "Employer" || group_df$to[1] == "Employer") {
        
        # Process each transition age range
        for(age_pair in key_ages) {
            # Find data points near transition
            nearby_data <- group_df %>%
                filter(age >= age_pair[1] - 1, age <= age_pair[2] + 1)
            
            if(nrow(nearby_data) > 0) {
                # Find indices in prediction grid near transition
                idx_range <- which(age_grid >= age_pair[1] & age_grid <= age_pair[2])
                
                if(length(idx_range) > 0) {
                    for(idx in idx_range) {
                        # Find nearest actual data point
                        nearest_idx <- which.min(abs(nearby_data$age - age_grid[idx]))
                        nearest_rate <- nearby_data$rate[nearest_idx]
                        
                        # Weight based on distance from transition midpoint
                        weight <- dnorm(age_grid[idx] - mean(age_pair), sd = 0.5)
                        predictions[idx] <- weight * nearest_rate + (1 - weight) * predictions[idx]
                    }
                }
            }
        }
    }
    
    return(list(
        predictions = predictions,
        age_grid = age_grid,
        df = group_df
    ))
}
# Update convert_to_tibble to handle the new output format
convert_to_tibble <- function(smoothed) {
    tibble(
        age = smoothed$age_grid,
        rate = if(!is.null(smoothed$predictions)) smoothed$predictions 
        else predict(smoothed$fit, data.frame(age = smoothed$age_grid)),
        from = first(smoothed$df$from),
        to = first(smoothed$df$to),
        model = first(smoothed$df$model),
        transition_type = first(smoothed$df$transition_type)
    )
}
rates_smoothed <- 
    grouped_data %>% map_df(~({
    .x %>% smooth_group_spline() %>% convert_to_tibble()
    }))

# Background mortality rates
mort_pre <- read_rds(here::here("results/hmd-mortality/hp_fit_2012.rds"))
mort_post <- read_rds(here::here("results/hmd-mortality/hp_fit_2019.rds"))

# aa = 22
# period = "pre"

construct_Q <- function(aa, period = "pre") {
    if (period=="pre") r_mort <- mort_pre(aa)
    if (period=="post") r_mort <- mort_post(aa)
    
    Q_ <- rates_smoothed %>% 
        mutate(rate = rate * 12) %>% 
        ungroup() %>% 
        filter(age ==aa & model == period) %>% 
        select(from,to,rate) %>% 
        spread(to,rate) %>% 
        column_to_rownames(var = "from") %>% 
        as.matrix()
    Q <- rbind(cbind(Q_, rep(r_mort, nrow(Q_))),rep(0,ncol(Q_)+1))
    dimnames(Q) = list(c(rownames(Q_),"Death"), c(colnames(Q_),"Death"))
    diag(Q) = -rowSums(Q,na.rm=TRUE)
    return(Q)
}


# Run a simple annual time step discrete time Markov model and 
# calcualte state occupancy by year. 

s0_pre <- 
    acs %>% 
    filter(age == params$age0) %>% 
    filter(time == "pre") %>% 
    pull(pct) 

s0_post <- 
    acs %>% 
    filter(age == params$age0) %>% 
    filter(time == "post") %>% 
    pull(pct)

plot_calibration <- function(rates_smoothed) {
    occ <- c(s0_pre,0)
    tr_pre <- params$ages_trace %>% map_df(~({
        Q_ <- construct_Q(.x, period = "pre")
        occ <<- occ %*% expm(Q_)
        data.frame(occ)
    })) %>%
        as_tibble() %>%
        mutate(surviving = Employer + OthPrivate + Public + Uninsured) %>%
        mutate(Employer = Employer / surviving,
               OthPrivate = OthPrivate / surviving,
               Public = Public / surviving,
               Uninsured = Uninsured / surviving) %>%
        select(-Death,-surviving) %>%
        mutate(time = "pre")  %>%
        mutate(age = params$ages_trace) %>%
        gather(type,pct,-time, -age)
    
    occ <- c(s0_post,0)
    tr_post <- params$ages_trace %>% map_df(~({
        Q_ <- construct_Q(.x, period = "post")
        occ <<- occ %*% expm(Q_)
        data.frame(occ)
    })) %>%
        as_tibble() %>%
        mutate(surviving = Employer + OthPrivate + Public + Uninsured) %>%
        mutate(Employer = Employer / surviving,
               OthPrivate = OthPrivate / surviving,
               Public = Public / surviving,
               Uninsured = Uninsured / surviving) %>%
        select(-Death,-surviving) %>%
        mutate(time = "post")  %>%
        mutate(age = params$ages_trace) %>%
        gather(type,pct,-time, -age)
    
    model <-
        tr_pre %>% bind_rows(tr_post)
    
    acs %>%
        ggplot(aes(x = age, y = pct))  + geom_point(colour = "red")  +
        ggthemes::theme_calc() + facet_wrap(time~type) +
        geom_point(data = model)
}

rates_smoothed %>% plot_calibration()


process_occupancy <- function(occ_df, period) {
    occ_df %>%
        as_tibble() %>%
        mutate(surviving = Employer + OthPrivate + Public + Uninsured) %>%
        mutate(
            Employer = Employer / surviving,
            OthPrivate = OthPrivate / surviving,
            Public = Public / surviving,
            Uninsured = Uninsured / surviving
        ) %>%
        select(-Death, -surviving) %>%
        mutate(
            time = period,
            age = params$ages_trace
        ) %>%
        gather(type, pct, -time, -age)
}


calibrate_rates_mcmc <- function(initial_rates, acs_targets, n_iterations = 1000, 
                                 burnin = 1000, thin = 10) {
    
    # Function to calculate model fit (negative log likelihood)
    calculate_fit <- function(rates) {
        rates_smoothed <<- rates  # Use global assignment for construct_Q
        
        # Run model for both periods
        occ <- c(s0_pre, 0)
        tr_pre <- params$ages_trace %>% map_df(~({
            Q_ <- construct_Q(.x, period = "pre")
            occ <<- occ %*% expm(Q_)
            data.frame(occ)
        })) %>% 
            process_occupancy("pre")
        
        occ <- c(s0_post, 0)
        tr_post <- params$ages_trace %>% map_df(~({
            Q_ <- construct_Q(.x, period = "post")
            occ <<- occ %*% expm(Q_)
            data.frame(occ)
        })) %>% 
            process_occupancy("post")
        
        model_pred <- bind_rows(tr_pre, tr_post)
        
        # Modified error calculation with weights
        total_error <- acs_targets %>%
            left_join(model_pred, by = c("age", "type", "time")) %>%
            mutate(
                # Add higher weight for OthPrivate
                weight = case_when(
                    type == "OthPrivate" ~ 2.0,  # Double weight for OthPrivate
                    TRUE ~ 1.0
                ),
                sq_error = (pct.x - pct.y)^2 * weight
            ) %>%
            ungroup() %>% 
            summarise(total_error = sum(sq_error, na.rm = TRUE)) %>%
            pull(total_error) %>%
            as.numeric() %>%
            `[`(1)
        
        return(total_error)
    }
    
    # Initialize
    current_rates <- initial_rates
    current_fit <- calculate_fit(current_rates)
    
    # Set up transition groups for targeted updates
    transition_groups <- current_rates %>%
        select(from, to, model) %>%
        distinct()
    
    # Store full chain history (after burnin, with thinning)
    chain_samples <- list()
    sample_counter <- 1
    
    # MCMC iterations
    for(iter in 1:n_iterations) {
        # Select a random transition group to update
        group_idx <- sample(1:nrow(transition_groups), 1)
        group <- transition_groups[group_idx,]
        
        # Create proposal
        proposed_rates <- current_rates
        mask <- proposed_rates$from == group$from & 
            proposed_rates$to == group$to & 
            proposed_rates$model == group$model
        
        scale_factor <- rnorm(1, 1, 0.1)
        proposed_rates$rate[mask] <- proposed_rates$rate[mask] * scale_factor
        
        # Calculate fit
        proposed_fit <- calculate_fit(proposed_rates)
        
        # Accept/reject
        accepted <- FALSE
        if(!is.na(proposed_fit) && (proposed_fit < current_fit || 
                                    runif(1) < exp(-(proposed_fit - current_fit)/0.01))) {
            current_rates <- proposed_rates
            current_fit <- proposed_fit
            accepted <- TRUE
        }
        
        # Store chain sample (after burnin and using thinning)
        if(iter > burnin && iter %% thin == 0) {
            chain_samples[[sample_counter]] <- list(
                rates = current_rates,
                fit = current_fit,
                iteration = iter
            )
            sample_counter <- sample_counter + 1
        }
        
        # Progress updates
        if(iter %% 100 == 0) {
            message("Completed iteration ", iter, " of ", n_iterations)
            message("Current fit: ", round(current_fit, 6))
            message("Samples collected: ", length(chain_samples))
        }
    }
    
    # Convert samples to a more useful format
    posterior_samples <- bind_rows(
        lapply(chain_samples, function(x) {
            x$rates %>% 
                mutate(
                    iteration = x$iteration,
                    fit = x$fit
                )
        })
    )
    
    # Calculate posterior summaries
    posterior_summaries <- posterior_samples %>%
        group_by(age, from, to, model) %>%
        summarise(
            rate_median = median(rate),
            rate_mean = mean(rate),
            rate_sd = sd(rate),
            rate_q025 = quantile(rate, 0.025),
            rate_q975 = quantile(rate, 0.975),
            .groups = 'drop'
        )
    
    # Return both point estimates and full posterior
    return(list(
        posterior_samples = posterior_samples,
        posterior_summaries = posterior_summaries,
        point_estimates = posterior_summaries %>% 
            select(age, from, to, model, rate = rate_median)
    ))
}

# Modified plotting function to show uncertainty
plot_calibration_with_uncertainty <- function(mcmc_results) {
    # Function to run model with specific rates
    run_model <- function(rates) {
        rates_smoothed <- rates
        
        # Pre period
        occ <- c(s0_pre, 0)
        tr_pre <- params$ages_trace %>% map_df(~({
            Q_ <- construct_Q(.x, period = "pre")
            occ <<- occ %*% expm(Q_)
            data.frame(occ)
        })) %>% process_occupancy("pre")
        
        # Post period
        occ <- c(s0_post, 0)
        tr_post <- params$ages_trace %>% map_df(~({
            Q_ <- construct_Q(.x, period = "post")
            occ <<- occ %*% expm(Q_)
            data.frame(occ)
        })) %>% process_occupancy("post")
        
        bind_rows(tr_pre, tr_post)
    }
    
    # Sample trajectories for uncertainty bands
    n_samples <- 100  # Number of trajectories to plot
    sample_indices <- sample(unique(mcmc_results$posterior_samples$iteration), n_samples)
    
    uncertainty_bands <- map_df(sample_indices, function(i) {
        rates <- mcmc_results$posterior_samples %>%
            filter(iteration == i) %>%
            select(-iteration, -fit)
        
        run_model(rates) %>%
            mutate(sample = i)
    })
    
    # Plot with uncertainty
    ggplot() +
        # Uncertainty bands
        stat_summary(data = uncertainty_bands,
                     aes(x = age, y = pct, group = sample),
                     geom = "line", alpha = 0.1) +
        # Median prediction
        geom_line(data = run_model(mcmc_results$point_estimates),
                  aes(x = age, y = pct)) +
        # Observed data
        geom_point(data = acs,
                   aes(x = age, y = pct),
                   color = "red") +
        facet_wrap(time ~ type) +
        ggthemes::theme_calc()
}

# Run MCMC with more iterations for posterior sampling
results <- calibrate_rates_mcmc(
    initial_rates = rates_smoothed,
    acs_targets = acs,
    n_iterations = 20000,  # More iterations for better posterior sampling
    burnin = 5000,        # Discard first 5000 iterations
    thin = 10             # Keep every 10th iteration
)

# Plot with uncertainty bands
plot_calibration_with_uncertainty(results)

# Save results
write_rds(results, here::here("results/calibrated-rates/calibrated-rates-with-uncertainty.rds"))
