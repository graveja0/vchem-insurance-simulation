# Prepare underlying MEPS data for multistate model
if (!exists("p_pre"))
    source(here::here("code/1_data_prep/download-and-prepare-meps-data.r"))
if (!exists("params"))
    source(here::here("code/0_setup/define-parameters.r"))
    
# Create skinny versions of pre and post datasets
create_skinny_df <- function(df) {
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
    
    # More granular age breaks based on observed transitions in the plot
    age_breaks <- c(
        0,                  # Start
        3,
        5,
        10,
        15,                 # Pre-transition
        18, 21,            # First major transition period
        23,                 # Post-college
        26, 28,            # Second major transition period
        35,                 # Early career
        45,                 # Mid career
        55,                 # Late career
        65                  # End
    )
    
    msm_data <- clean_data %>%
        mutate(
            state = state_to_num[type],
            age_group = cut(age, breaks = age_breaks, right = FALSE)
        )
    
    # Store attributes
    attr(msm_data, "states") <- states
    attr(msm_data, "state_to_num") <- state_to_num
    attr(msm_data, "age_breaks") <- age_breaks
    
    return(msm_data)
}

# Create skinny versions of both datasets
df_skinny_pre <- df <- create_skinny_df(df_final_pre)
df_skinny_post <- create_skinny_df(df_final_post)

fit_multistate_model <- function(df) {
    # Find weight column
    weight_col <- names(df)[tolower(names(df)) == "longwt"]
    if(length(weight_col) == 0) {
        stop("Weight column not found. Please ensure 'longwt' or 'LONGWT' exists in the data.")
    }
    
    all_attrs <- attributes(df)
    age_breaks <- attr(df, "age_breaks")
    
    # Print weighted transition patterns
    transition_summary <- df %>%
        group_by(age_group) %>%
        arrange(dupersid, time) %>%
        group_by(age_group, dupersid) %>%
        summarise(
            transitions = paste(state, collapse = "->"),
            weight = mean(!!sym(weight_col)),
            .groups = "drop"
        ) %>%
        group_by(age_group, transitions) %>%
        summarise(
            n = n(),
            weighted_n = sum(weight),
            .groups = "drop"
        ) %>%
        arrange(age_group, desc(weighted_n))
    
    message("Weighted transition patterns by age group:")
    print(transition_summary)
    
    # Print weighted sample sizes
    sample_sizes <- df %>%
        group_by(age_group) %>%
        summarise(
            n_people = n_distinct(dupersid),
            weighted_n_people = sum(!!sym(weight_col)[!duplicated(dupersid)]),
            n_obs = n(),
            weighted_n_obs = sum(!!sym(weight_col)),
            n_transitions = n() - n_distinct(dupersid),
            weighted_n_transitions = sum(!!sym(weight_col)) - 
                sum(!!sym(weight_col)[!duplicated(dupersid)])
        )
    
    message("\nWeighted sample sizes by age group:")
    print(sample_sizes)
    
    # Fit separate models for each age group
    models <- list()
    
    for(ag in unique(df$age_group)) {
        message("\nFitting model for age group: ", ag)
        
        df_sub <- df %>% 
            filter(age_group == ag) %>%
            group_by(dupersid) %>%
            filter(n() > 1) %>%
            ungroup()
        
        if(nrow(df_sub) > 0) {
            # Calculate initial Q matrix using weighted transitions
            Q_init <- calculate_Q_init(df_sub, weight_var = weight_col)
            
            models[[as.character(ag)]] <- tryCatch({
                msm(
                    state ~ time,
                    subject = dupersid,
                    data = df_sub,
                    qmatrix = Q_init,
                    method = "BFGS",
                    control = list(
                        fnscale = 1000,
                        maxit = 1000,
                        reltol = 1e-4,
                        trace = 1,
                        REPORT = 1
                    )
                )
            }, error = function(e) {
                message("Error fitting model for age group ", ag, ": ", e$message)
                return(NULL)
            })
        }
    }
    
    # Store attributes and weighted summaries
    for (attr_name in names(all_attrs)) {
        if (!attr_name %in% c("names", "class", "row.names")) {
            attr(models, attr_name) <- all_attrs[[attr_name]]
        }
    }
    attr(models, "sample_sizes") <- sample_sizes
    attr(models, "transition_summary") <- transition_summary
    attr(models, "weighted_init") <- TRUE  # Flag that we used weighted initialization
    
    return(models)
}

# Update Q matrix initialization to use dynamic weight column
calculate_Q_init <- function(df, weight_var) {
    n_states <- length(unique(df$state))
    
    # Calculate weighted transition counts and total time
    transitions <- df %>%
        group_by(dupersid) %>%
        arrange(time) %>%
        summarise(
            from = head(state, -1),
            to = tail(state, -1),
            weight = mean(!!sym(weight_var)),  # Use dynamic column name
            time_diff = diff(time)
        ) %>%
        ungroup()
    
    # Calculate weighted total time spent in each state
    total_time <- df %>%
        group_by(state) %>%
        summarise(total_time = sum(!!sym(weight_var))) %>%  # Use dynamic column name
        pull(total_time)
    
    # Create transition matrix
    Q <- matrix(0, n_states, n_states)
    
    # Calculate weighted rates
    for(i in 1:n_states) {
        for(j in 1:n_states) {
            if(i != j) {
                weighted_transitions <- sum(transitions$weight[transitions$from == i & transitions$to == j])
                if(total_time[i] > 0) {
                    Q[i,j] <- weighted_transitions / total_time[i]
                }
            }
        }
    }
    
    # Ensure diagonal elements make rows sum to 0
    diag(Q) <- -rowSums(Q)
    
    # Add small positive value to prevent numerical issues
    Q[Q == 0] <- 1e-6
    diag(Q) <- -rowSums(Q)
    
    return(Q)
}



