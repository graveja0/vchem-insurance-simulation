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



