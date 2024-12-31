# Prepare underlying MEPS data for multistate model
if (!exists("p_pre"))
    source(here::here("code/1_data_prep/download-and-prepare-meps-data.r"))
if (!exists("params"))
    source(here::here("code/0_setup/define-parameters.r"))
    
knot_locations <- c(0,1,5, 10, 15, 18, 19, 20, 21, 26, 35, 50, 60)
knot_locations <- knot_locations[between(knot_locations,params_$age0,params_$age0+params_$horizon)]

# Create skinny versions of pre and post datasets
create_skinny_df <- function(df, 
                           # Place more basis functions around key transition ages
                           centers = c(
                               seq(0, 15, by = 5),    # Childhood changes
                               seq(16, 22, by = 1),   # Dense coverage of 18-20 transition
                               seq(23, 29, by = 1),   # Dense coverage of 26 transition
                               seq(30, 60, by = 5)    # Sparser coverage of stable periods
                           ), 
                           width = 2) {  # Narrower width to capture sharp transitions
    
    # Keep the original cleaning code exactly as specified
    clean_data <- df %>%
        filter(age >= params$age0 & age <= (params$age0 + params$horizon)) %>% 
        group_by(dupersid) %>%
        filter(!is.na(type)) %>% 
        filter(n() > 2) %>% 
        arrange(dupersid, age, desc(row_number())) %>%
        distinct(dupersid, month, .keep_all = TRUE) %>%
        arrange(dupersid, month) %>%
        group_by(dupersid) %>%
        mutate(time = row_number()) %>%
        ungroup()
    
    # Ensure states are numeric for msm
    states <- sort(unique(clean_data$type))
    state_to_num <- seq_along(states)
    names(state_to_num) <- states
    
    # Create adaptive width Gaussian basis functions
    gaussian_terms <- sapply(centers, function(c) {
        # Use smaller width for transition ages
        local_width <- if(c >= 16 && c <= 29) width else width * 2
        exp(-0.5 * ((clean_data$age - c)/local_width)^2)
    })
    colnames(gaussian_terms) <- paste0("age_gauss", centers)
    
    # Add polynomial terms for long-term trends
    poly_terms <- poly(clean_data$age, degree = 2, raw = TRUE)
    colnames(poly_terms) <- paste0("age_poly", 1:2)
    
    # Combine Gaussian and polynomial terms
    combined_terms <- cbind(gaussian_terms, poly_terms)
    
    # Scale all terms to [0,1]
    scaled_terms <- apply(combined_terms, 2, function(x) {
        (x - min(x)) / (max(x) - min(x))
    })
    
    basis_data <- as.data.frame(scaled_terms)
    
    # Combine with original data
    msm_data <- clean_data %>%
        mutate(
            state = state_to_num[type]
        ) %>%
        bind_cols(basis_data)
    
    # Store attributes
    attr(msm_data, "states") <- states
    attr(msm_data, "state_to_num") <- state_to_num
    attr(msm_data, "centers") <- centers
    attr(msm_data, "width") <- width
    attr(msm_data, "basis_range") <- apply(combined_terms, 2, range)
    
    return(msm_data)
}

# Function to create prediction data for a single age
create_prediction_data <- function(age, model) {
    centers <- attr(model, "centers")
    width <- attr(model, "width")
    basis_range <- attr(model, "basis_range")
    
    # Create Gaussian basis values with adaptive width
    gaussian_terms <- sapply(centers, function(c) {
        local_width <- if(c >= 16 && c <= 29) width else width * 2
        exp(-0.5 * ((age - c)/local_width)^2)
    })
    
    # Add polynomial terms
    poly_terms <- poly(age, degree = 2, raw = TRUE)
    combined_terms <- c(gaussian_terms, poly_terms)
    
    # Scale using stored ranges
    scaled_terms <- mapply(function(x, ranges) {
        (x - ranges[1]) / (ranges[2] - ranges[1])
    }, combined_terms, split(basis_range, col(basis_range)))
    
    pred_data <- as.list(scaled_terms)
    names(pred_data) <- c(paste0("age_gauss", centers), paste0("age_poly", 1:2))
    
    return(pred_data)
}

fit_multistate_model <- function(df, method = "BFGS") {
    all_attrs <- attributes(df)
    
    # Combine Gaussian and polynomial terms
    basis_cols <- c(
        grep("^age_gauss", names(df), value = TRUE),
        grep("^age_poly", names(df), value = TRUE)
    )
    covariate_formula <- paste(basis_cols, collapse = " + ")
    
    # Very conservative initial values
    n_states <- length(attr(df, "states"))
    Q <- matrix(1e-8, n_states, n_states)
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
    
    # Fit model with conservative settings
    msm_fit <- msm(
        state ~ time,
        subject = dupersid,
        data = df,
        qmatrix = Q,
        covariates = formula(paste("~", covariate_formula)),
        control = list(
            fnscale = 50000,
            maxit = 1000,
            reltol = 1e-5,
            trace = 1,
            REPORT = 1
        ),
        method = method,
        center = TRUE
    )
    
    # Store attributes
    for (attr_name in names(all_attrs)) {
        if (!attr_name %in% c("names", "class", "row.names")) {
            attr(msm_fit, attr_name) <- all_attrs[[attr_name]]
        }
    }
    
    return(msm_fit)
}

# Function to get transition intensities for a specific age
get_Q_matrix <- function(msm_fit, age) {
    # Create prediction data
    pred_data <- create_prediction_data(age, msm_fit)
    
    # Get Q matrix
    states <- attr(msm_fit, "states")
    n_states <- length(states)
    
    # Get raw qmatrix values
    q_values <- qmatrix.msm(msm_fit, covariates = pred_data)
    
    # Convert to matrix form
    Q <- matrix(0, n_states, n_states)
    for(i in 1:n_states) {
        for(j in 1:n_states) {
            if(i != j) {
                Q[i,j] <- as.numeric(q_values[i,j][1])
            }
        }
    }
    
    diag(Q) <- -rowSums(Q)
    rownames(Q) <- colnames(Q) <- states
    
    return(Q)
}

# Create skinny versions of both datasets
df_skinny_pre <- df <- create_skinny_df(df_final_pre, knot_locations)
df_skinny_post <- create_skinny_df(df_final_post, knot_locations)

# Fit models on pre and post datasets
msm_fit_pre <- fit_multistate_model(df_skinny_pre)
# msm_fit_post <- fit_multistate_model(df_skinny_post)

# Example: Get Q matrix for age 22
Q_22 <- get_Q_matrix(msm_fit_pre, age = 22)
print("Transition intensity matrix for age 22:")
print(round(expm(Q_22), 4))

