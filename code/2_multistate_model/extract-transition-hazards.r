if (!exists("p_pre"))
    source(here::here("code/1_data_prep/download-and-prepare-meps-data.r"))
if (!exists("params"))
    source(here::here("code/0_setup/define-parameters.r"))
if (!exists("df_skinny_post"))
    source(here::here("code/1_data_prep/prepare-meps-data-for-multistate-model.r"))
if (!exists("msm_fit_pre"))
    source(here::here("code/2_multistate_model/fit-multistate-model.r"))


# Function to extract transition intensity matrix for a given age and model fit
get_Q_matrix <- function(models, age) {
    age_breaks <- attr(models, "age_breaks")
    
    # Find appropriate age group
    for(i in 1:(length(age_breaks)-1)) {
        if(age >= age_breaks[i] && age < age_breaks[i+1]) {
            ag <- cut(age, breaks = age_breaks, right = FALSE)
            model <- models[[as.character(ag)]]
            if(!is.null(model)) {
                Q <- qmatrix.msm(model,ci="none")
                # Convert to matrix format
                Q_mat <- matrix(0, nrow=nrow(Q), ncol=ncol(Q))
                Q_mat[] <- Q
                return(Q_mat)
            }
        }
    }
    return(NULL)
}

# Add this function after the existing functions

process_msm_results <- function(msm_fit_pre, msm_fit_post, params) {
    # Create sequence of ages to evaluate
    ages <- params$age0:min(65, params$age0 + params$horizon)
    
    # Get state mappings from attributes
    states <- attr(msm_fit_pre, "states")
    state_to_num <- attr(msm_fit_pre, "state_to_num")
    num_to_state <- names(state_to_num)
    names(num_to_state) <- state_to_num
    
    # Function to extract transition rates from Q matrix
    extract_rates <- function(Q) {
        n_states <- nrow(Q)
        transitions <- expand.grid(from_num = 1:n_states, to_num = 1:n_states) %>%
            filter(from_num != to_num) %>%
            mutate(
                from = factor(num_to_state[as.character(from_num)], 
                             levels = seq_along(params$v_tr_names),
                             labels = params$v_tr_names),
                to = factor(num_to_state[as.character(to_num)], 
                          levels = seq_along(params$v_tr_names),
                          labels = params$v_tr_names),
                transition_type = paste0(from_num, "->", to_num),
                rate = sapply(seq_len(nrow(.)), function(i) Q[from_num[i], to_num[i]])
            ) %>%
            select(from, to, transition_type, rate)
        return(transitions)
    }
    
    # Process each age for both models
    results <- data.frame()
    
    for(age in ages) {
        # Pre-ACA rates
        Q_pre <- get_Q_matrix(msm_fit_pre, age)
        if(!is.null(Q_pre)) {
            pre_rates <- extract_rates(Q_pre) %>%
                mutate(
                    age = age,
                    model = "pre"
                )
            results <- bind_rows(results, pre_rates)
        }
        
        # Post-ACA rates
        Q_post <- get_Q_matrix(msm_fit_post, age)
        if(!is.null(Q_post)) {
            post_rates <- extract_rates(Q_post) %>%
                mutate(
                    age = age,
                    model = "post"
                )
            results <- bind_rows(results, post_rates)
        }
    }
    
    # Ensure consistent column order
    results <- results %>%
        select(age, from, to, transition_type, model, rate) %>%
        arrange(model, age, from, to)
    
    return(results)
}



# Example usage:
transition_rates <- process_msm_results(msm_fit_pre, msm_fit_post, params) %>% 
    filter(!is.na(to) & !is.na(from))
saveRDS(transition_rates, here::here("results/model_objects/transition_rates.rds"))


