# Prepare underlying MEPS data for multistate model
if (!exists("p_pre"))
    source(here::here("code/1_data_prep/download-and-prepare-meps-data.r"))
if (!exists("params"))
    source(here::here("code/0_setup/define-parameters.r"))
    
knot_locations <- c(0,1,5, 10, 15, 18, 19, 20, 21, 26, 35, 50, 60)
knot_locations <- knot_locations[between(knot_locations,params_$age0,params_$age0+params_$horizon)]

# Create skinny versions of pre and post datasets
create_skinny_df <- function(df, knots) {
    # Clean data first - more aggressive cleaning
    clean_data <- df %>%
        filter(age >= params$age0 & age <= params_$age0+params_$horizon) %>% 
        filter(!is.na(type)) %>% 
        # Remove subjects with only one observation
        group_by(dupersid) %>%
        arrange(dupersid, month) %>% 
        mutate(time = row_number()) %>% 
        filter(n() > 2) %>%  # Require at least 3 observations
        distinct(dupersid, time, .keep_all = TRUE) %>%
        # Ensure observations are properly ordered
        arrange(dupersid, time) 
        
    
    # Ensure states are numeric for msm
    states <- sort(unique(clean_data$type))
    state_to_num <- seq_along(states)
    names(state_to_num) <- states
    
    # Create spline basis with specified knots
    ns_basis <- ns(clean_data$age, 
                   knots = knots[-c(1, length(knots))],
                   Boundary.knots = c(knots[1], knots[length(knots)]))
    
    # Get spline basis values and scale them
    spline_matrix <- predict(ns_basis, clean_data$age)
    colnames(spline_matrix) <- paste0("age_spline", 1:ncol(spline_matrix))
    
    # Store the mapping of age to spline values
    age_to_splines <- clean_data %>%
        select(age) %>%
        bind_cols(as.data.frame(spline_matrix)) %>%
        distinct(age, .keep_all = TRUE) %>%
        arrange(age)
    
    # Create the model matrix for splines
    spline_data <- as.data.frame(spline_matrix)
    
    # Combine with original data
    msm_data <- clean_data %>%
        mutate(
            state = state_to_num[type]  # Convert states to numeric
        ) %>%
        bind_cols(spline_data)
    
    # Create function to get scaled splines for a given age
    get_scaled_splines <- function(age) {
        # Get raw spline values
        raw_splines <- predict(ns_basis, age)
        # Scale using same parameters as original data
        scaled_splines <- scale(raw_splines, 
                              center = attr(spline_matrix, "scaled:center"),
                              scale = attr(spline_matrix, "scaled:scale"))
        # Convert to vector and name elements
        spline_vec <- as.vector(scaled_splines)
        names(spline_vec) <- paste0("age_spline", 1:length(spline_vec))
        return(spline_vec)
    }
    
    # Store attributes
    attr(msm_data, "states") <- states
    attr(msm_data, "state_to_num") <- state_to_num
    attr(msm_data, "ns_basis") <- ns_basis
    attr(msm_data, "knots") <- knots
    attr(msm_data, "spline_names") <- colnames(spline_matrix)
    attr(msm_data, "age_to_splines") <- age_to_splines
    attr(msm_data, "spline_scaling") <- list(
        center = attr(spline_matrix, "scaled:center"),
        scale = attr(spline_matrix, "scaled:scale")
    )
    attr(msm_data, "get_scaled_splines") <- get_scaled_splines
    
    return(msm_data)
}

# Create skinny versions of both datasets
df_skinny_pre <- df <- create_skinny_df(df_final_pre, knot_locations)
df_skinny_post <- create_skinny_df(df_final_post, knot_locations)

# Function to fit multistate model using msm package
fit_multistate_model <- function(df, method = "BFGS") {
    # Get all attributes from the input data
    all_attrs <- attributes(df)
    
    # Create formula with all spline terms
    spline_cols <- grep("^age_spline", names(df), value = TRUE)
    covariate_formula <- paste(spline_cols, collapse = " + ")
    
    # Create transition matrix for initial values with VERY small values
    n_states <- length(attr(df, "states"))
    Q <- matrix(0.0000001, n_states, n_states)  # Much smaller initial values
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
    
    # First try fitting with just age
    msm_fit <- msm(
        state ~ time,
        subject = dupersid,
        data = df,
        qmatrix = Q,
        covariates = ~age,  # Start with just age
        control = list(
            fnscale = 10000,     # Increased scaling
            maxit = 1000,
            reltol = 1e-8,
            trace = 1,
            REPORT = 1
        ),
        method = method,  # Use specified method
        center = TRUE
    )
    
    # Use the converged Q matrix as initial values for the spline model
    Q_init <- qmatrix.msm(msm_fit)
    
    # Now fit the full spline model
    msm_fit <- msm(
        state ~ time,
        subject = dupersid,
        data = df,
        qmatrix = Q_init,
        covariates = formula(paste("~", covariate_formula)),
        control = list(
            fnscale = 10000,
            maxit = 1000,
            reltol = 1e-8,
            trace = 1,
            REPORT = 1
        ),
        method = method,  # Use specified method
        center = TRUE
    )
    
    # Store all attributes from the original data
    for (attr_name in names(all_attrs)) {
        if (!attr_name %in% c("names", "class", "row.names")) {
            attr(msm_fit, attr_name) <- all_attrs[[attr_name]]
        }
    }
    
    # Explicitly store key attributes
    msm_fit$data <- df
    msm_fit$states <- attr(df, "states")
    msm_fit$state_to_num <- attr(df, "state_to_num")
    msm_fit$knots <- attr(df, "knots")
    msm_fit$spline_names <- attr(df, "spline_names")
    msm_fit$ns_basis <- attr(df, "ns_basis")
    msm_fit$age_to_splines <- attr(df, "age_to_splines")
    msm_fit$spline_scaling <- attr(df, "spline_scaling")
    
    return(msm_fit)
}

# Fit models on pre and post datasets
msm_fit_pre <- fit_multistate_model(df_skinny_pre)
# msm_fit_post <- fit_multistate_model(df_skinny_post)

# Function to extract transition intensity matrix for a given age and model fit
get_Q_matrix <- function(msm_fit, age, type = "predicted") {
    # Get the age-to-splines mapping
    age_to_splines <- attr(msm_fit$data, "age_to_splines")
    
    # Find the exact spline values used for this age in the training data
    spline_values <- age_to_splines %>%
        filter(age == !!age) %>%
        select(starts_with("age_spline")) %>%
        as.list()
    
    # If age not found in original data, compute new spline values
    if (length(spline_values) == 0) {
        ns_basis <- attr(msm_fit$data, "ns_basis")
        raw_splines <- predict(ns_basis, age)
        scaling <- attr(msm_fit$data, "spline_scaling")
        scaled_splines <- scale(matrix(raw_splines, ncol = length(raw_splines)), 
                              center = scaling$center,
                              scale = scaling$scale)[1,]
        spline_values <- as.list(setNames(scaled_splines, 
                                        paste0("age_spline", seq_along(scaled_splines))))
    }
    
    # Get Q matrix
    states <- attr(msm_fit, "states")
    n_states <- length(states)
    
    # Get raw qmatrix values
    q_values <- qmatrix.msm(msm_fit)
    
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

# Example: Get Q matrix for age 22
Q_22 <- get_Q_matrix(msm_fit_pre, age = 22)
print("Transition intensity matrix for age 22:")
print(round(expm(Q_22), 4))

