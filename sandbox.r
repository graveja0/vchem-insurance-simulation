library(tidyverse)
library(msm)

# Function to fit piecewise models with more reasonable age groups
fit_piecewise_model <- function(df, 
                              # Wider age groups, especially around transitions
                              age_breaks = c(0, 16, 18, 22, 25, 30, 40, 50, 65)) {
    
    # Add age group variable
    df <- df %>% 
        mutate(age_group = cut(age, breaks = age_breaks, right = FALSE))
    
    # Check sample sizes before fitting
    sample_sizes <- df %>%
        group_by(age_group) %>%
        summarise(
            n_people = n_distinct(dupersid),
            n_obs = n(),
            n_transitions = n() - n_distinct(dupersid)
        )
    
    message("Sample sizes by age group:")
    print(sample_sizes)
    
    # Fit separate models for each age group
    models <- list()
    for(ag in unique(df$age_group)) {
        message("\nFitting model for age group: ", ag)
        
        df_sub <- df %>% 
            filter(age_group == ag) %>%
            # Ensure we have enough transitions in this group
            group_by(dupersid) %>%
            filter(n() > 1) %>%
            ungroup()
        
        if(nrow(df_sub) > 0) {
            # Conservative initial values
            Q <- matrix(0.0001, 4, 4)
            diag(Q) <- 0
            diag(Q) <- -rowSums(Q)
            
            models[[as.character(ag)]] <- tryCatch({
                msm(
                    state ~ time,
                    subject = dupersid,
                    data = df_sub,
                    qmatrix = Q,
                    control = list(
                        fnscale = 4000,
                        maxit = 1000,
                        reltol = 1e-6
                    )
                )
            }, error = function(e) {
                message("Error fitting model for age group ", ag, ": ", e$message)
                return(NULL)
            })
        }
    }
    
    # Add metadata
    attr(models, "age_breaks") <- age_breaks
    attr(models, "sample_sizes") <- sample_sizes
    return(models)
}

# Function to get Q matrix for any age
get_Q_matrix <- function(models, age) {
    age_breaks <- attr(models, "age_breaks")
    
    # Find appropriate age group
    for(i in 1:(length(age_breaks)-1)) {
        if(age >= age_breaks[i] && age < age_breaks[i+1]) {
            ag <- cut(age, breaks = age_breaks, right = FALSE)
            model <- models[[as.character(ag)]]
            if(!is.null(model)) {
                return(qmatrix.msm(model))
            }
        }
    }
    return(NULL)
}

# Test the approach
df <- readRDS("data/meps_cleaned.rds")
piecewise_models <- fit_piecewise_model(df)

# Function to examine results
examine_transitions <- function(models, ages = c(17, 20, 24, 27, 35, 45)) {
    results <- list()
    for(age in ages) {
        Q <- get_Q_matrix(models, age)
        if(!is.null(Q)) {
            results[[as.character(age)]] <- Q
        }
    }
    return(results)
}

# Plot transition probabilities over time for each initial state
plot_transition_probs <- function(models, from_state = 1, t = 1) {
    age_breaks <- attr(models, "age_breaks")
    ages <- seq(min(age_breaks), max(age_breaks)-1, by = 1)
    
    probs <- data.frame()
    
    for(age in ages) {
        Q <- get_Q_matrix(models, age)
        if(!is.null(Q)) {
            P <- MatrixExp(Q * t)  # 1-year transition probabilities
            probs <- rbind(probs,
                          data.frame(
                              age = age,
                              to_state = 1:4,
                              probability = P[from_state,]
                          ))
        }
    }
    
    ggplot(probs, aes(x = age, y = probability, color = factor(to_state))) +
        geom_line() +
        theme_minimal() +
        labs(title = paste("1-year transition probabilities from state", from_state),
             color = "To State")
}

# Save results
results <- examine_transitions(piecewise_models)
saveRDS(piecewise_models, "output/piecewise_models.rds")
saveRDS(results, "output/transition_matrices.rds")

# Generate plots for each initial state
plots <- lapply(1:4, function(s) plot_transition_probs(piecewise_models, from_state = s)) 