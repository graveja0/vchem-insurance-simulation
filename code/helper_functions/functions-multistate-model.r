




prepare_msm_data <- function(data, knot_locations) {
    # Ensure states are numeric for msm
    states <- sort(unique(data$type))
    state_to_num <- seq_along(states)
    names(state_to_num) <- states
    
    # Create spline basis with specified knots
    age_range <- range(data$age)
    ns_basis <- ns(data$age, 
                   knots = knot_locations[-c(1, length(knot_locations))],
                   Boundary.knots = c(knot_locations[1], knot_locations[length(knot_locations)]))
    
    # Get spline basis values
    spline_matrix <- predict(ns_basis, data$age)
    colnames(spline_matrix) <- paste0("age_spline", 1:ncol(spline_matrix))
    
    # Create the model matrix for splines and scaled weight
    spline_data <- as.data.frame(spline_matrix)
    
    # Scale weight to mean 0, sd 1
    weight_scaled <- scale(data$weight)
    spline_data$weight_scaled <- as.vector(weight_scaled)
    
    # Combine with original data
    msm_data <- data %>%
        mutate(
            state = state_to_num[type]  # Convert states to numeric
        ) %>%
        bind_cols(spline_data)
    
    return(list(
        data = msm_data,
        states = states,
        state_to_num = state_to_num,
        ns_basis = ns_basis,
        knots = knot_locations,
        spline_names = c(colnames(spline_matrix), "weight_scaled"),
        weight_mean = mean(data$weight),
        weight_sd = sd(data$weight)
    ))
}


#' Fit multi-state model with spline-based age effects
fit_insurance_msm <- function(prepared_data, control = list(maxit = 1000, reltol = 1e-8, fnscale = 4000)) {
    # Create formula for spline terms and weight
    covariates <- prepared_data$spline_names  # Now includes weight_scaled
    covariate_formula <- paste(covariates, collapse = " + ")
    
    # Create transition matrix for initial values
    n_states <- length(prepared_data$states)
    Q <- matrix(0.1, n_states, n_states)
    diag(Q) <- -rowSums(Q)
    
    # Create a simple constraint structure
    # Format: list of numeric vectors
    n_params <- length(covariates)
    
    constraint <- list(
        # Constraint for Public (3) to Employer (1) transitions
        c(rep(1, n_params)),  # All coefficients positive
        
        # Constraint for Public (3) to OthPrivate (2) transitions
        c(rep(1, n_params)),  # All coefficients positive
        
        # Constraint for Employer (1) to Public (3) transitions
        c(rep(-1, n_params))  # All coefficients negative
    )
    
    # Fit model
    model <- msm(
        state ~ month,
        subject = patient_id,
        data = prepared_data$data,
        qmatrix = Q,
        covariates = formula(paste("~", covariate_formula)),
        control = control,
        method = "BFGS",
        constraint = constraint
    )
    
    # Add useful attributes
    attr(model, "states") <- prepared_data$states
    attr(model, "state_to_num") <- prepared_data$state_to_num
    attr(model, "knots") <- prepared_data$knots
    attr(model, "spline_names") <- prepared_data$spline_names
    
    return(model)
}
