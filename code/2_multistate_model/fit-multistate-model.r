if (!exists("p_pre"))
    source(here::here("code/1_data_prep/download-and-prepare-meps-data.r"))
if (!exists("params"))
    source(here::here("code/0_setup/define-parameters.r"))
if (!exists("df_skinny_post"))
    source(here::here("code/1_data_prep/prepare-meps-data-for-multistate-model.r"))

# Helper function to calculate empirical transition rates
calculate_Q_init <- function(df) {
    n_states <- length(unique(df$state))
    
    # Calculate transition counts and total time in each state
    transitions <- df %>%
        group_by(dupersid) %>%
        arrange(time) %>%
        summarise(
            from = head(state, -1),
            to = tail(state, -1),
            time_diff = diff(time)
        ) %>%
        ungroup()
    
    # Calculate total time spent in each state
    total_time <- df %>%
        group_by(state) %>%
        summarise(total_time = n()) %>%
        pull(total_time)
    
    # Create transition matrix
    Q <- matrix(0, n_states, n_states)
    
    # Calculate rates
    for(i in 1:n_states) {
        for(j in 1:n_states) {
            if(i != j) {
                n_transitions <- sum(transitions$from == i & transitions$to == j)
                if(total_time[i] > 0) {
                    Q[i,j] <- n_transitions / total_time[i]
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

# Function to fit multistate model using msm package
fit_multistate_model <- function(df) {
    all_attrs <- attributes(df)
    age_breaks <- attr(df, "age_breaks")
    
    # Analyze transitions
    transition_summary <- df %>%
        group_by(age_group) %>%
        arrange(dupersid, time) %>%
        group_by(age_group, dupersid) %>%
        summarise(
            transitions = paste(state, collapse = "->"),
            .groups = "drop"
        ) %>%
        count(age_group, transitions) %>%
        arrange(age_group, desc(n))
    
    message("Transition patterns by age group:")
    print(transition_summary)
    
    # Print sample sizes
    sample_sizes <- df %>%
        group_by(age_group) %>%
        summarise(
            n_people = n_distinct(dupersid),
            n_obs = n(),
            n_transitions = n() - n_distinct(dupersid)
        )
    
    message("\nSample sizes by age group:")
    print(sample_sizes)
    
    # Fit models
    models <- list()
    
    for(ag in unique(df$age_group)) {
        message("\nFitting model for age group: ", ag)
        
        df_sub <- df %>% 
            filter(age_group == ag) %>%
            group_by(dupersid) %>%
            filter(n() > 1) %>%
            ungroup()
        
        if(nrow(df_sub) > 0) {
            # Calculate initial Q matrix from observed transitions
            Q_init <- calculate_Q_init(df_sub)
            
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
    
    # Store attributes
    for (attr_name in names(all_attrs)) {
        if (!attr_name %in% c("names", "class", "row.names")) {
            attr(models, attr_name) <- all_attrs[[attr_name]]
        }
    }
    attr(models, "sample_sizes") <- sample_sizes
    attr(models, "transition_summary") <- transition_summary
    
    return(models)
}

# Fit models on pre and post datasets
msm_fit_pre <- fit_multistate_model(df_skinny_pre)
msm_fit_post <- fit_multistate_model(df_skinny_post)


