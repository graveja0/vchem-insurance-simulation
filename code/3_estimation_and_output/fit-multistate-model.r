if (!exists("p_pre"))
    source(here::here("R/download-and-prepare-meps-data.r"))
if (!exists("params"))
    source(here::here("R/define-parameters.r"))
source(here::here("R/functions-multistate-model.r"))

knot_locations <- c(0,1,5, 10, 15, 18, 19, 20, 21, 26, 35, 50, 60)
knot_locations <- knot_locations[between(knot_locations,params_$age0,params_$age0+params_$horizon)]

####################
# ACS-Based Targets
####################

targets <- 
    read_rds(here::here(glue::glue("_inputs/acs-calibration-targets/acs-calibration-targets_2012.rds"))) %>% 
    mutate(type = case_when(type==1 ~ "Employer",
                            type==2 ~ "OthPrivate",
                            type==3 ~ "Public", 
                            type==4 ~ "Uninsured")) %>% 
    filter(agep >= params$age0 & agep <= params$age0 + params$horizon)

targets %>% 
    ggplot(aes(x = agep, y = pct)) + geom_line() + facet_wrap(~type) + 
    ggthemes::theme_economist_white()


targets %>% 
    ggplot(aes(x = agep, y = pct)) + 
    geom_line() + 
    geom_vline(xintercept = knot_locations, 
               linetype = "dashed", 
               color = "red", 
               alpha = 0.5) +
    facet_wrap(~type) + 
    ggthemes::theme_economist_white() +
    labs(title = "Insurance Coverage by Age",
         subtitle = "Red dashed lines show spline knot locations",
         x = "Age",
         y = "Percentage")

####################
# MEPS Input Data
####################

df =
    df_final_pre %>%
    group_by(dupersid) %>%
    mutate(month = (year - min(as.numeric(pre_lut))) * 12 + month(month)) %>%
    # Exclusion: no nonelderly medicare
    mutate(nonelderly_medicare = max(type == 5 & age < 65)) %>%
    filter(nonelderly_medicare == 0) %>%
    group_by(dupersid) %>%
    # Exclusion: 3+ months observed
    filter(n() >= 3)  %>%
    # Exclusion: age restriciton.
    filter(age >= params$age0 & age <= params$age0 + params$horizon) %>% 
    mutate(type = paste0(factor(
        type, levels = 1:4, labels = params$v_tr_names
    ))) %>%
    select(patient_id = dupersid,
           weight = longwt,
           month,
           type,
           age = age)

set.seed(123)
M = 1000
sampled_ids <- sample(unique(df$patient_id), M)
df_model <- df %>% filter(patient_id %in% sampled_ids)

# Fit model
# First fit model with more knots
model <- fit_insurance_msm(
    prepared_data = prepare_msm_data(
        data = df_model,
        knot_locations = knot_locations  # More detailed knots
    )
)

get_qmatrix <- function(age, model, prepared_data) {
    # Get spline basis values for this age
    spline_vals <- predict(prepared_data$ns_basis, age)
    
    # Create data frame with spline values and mean weight
    pred_data <- data.frame(matrix(spline_vals, nrow=1))
    colnames(pred_data) <- paste0("age_spline", 1:ncol(spline_vals))
    pred_data$weight_scaled <- 0  # Use mean weight (since scaled, mean is 0)
    
    # Get Q matrix from msm package
    qmatrix <- msm::qmatrix.msm(model, covariates = pred_data)
    
    return(qmatrix)
}

get_qmatrix(
    age = 18,
    model = model,
    prepared_data = prepare_msm_data(data = df_model, knot_locations = knot_locations)  # More detailed knots
)

    