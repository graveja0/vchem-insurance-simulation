if (!exists("manifest")) 
    source(here::here("code/0_setup/manifest.r"))
if (!exists("params"))
    source(here::here("code/0_setup/define-parameters.r"))

# ACS Targets

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


transition_rates <-  # created in examine-model-fit.r
    read_rds(here::here("results/model_objects/transition_rates.rds"))

grouped_data <- 
    transition_rates %>%
    group_by(from, to, model) %>%
    group_split()

k = 20

# Function to smooth one group
smooth_group_spline <- function(group_df) {
    # Create prediction grid
    age_grid <- seq(min(group_df$age), max(group_df$age), by = 0.5)
    
    # Fit GAM with cubic splines
    gam_fit <- gam(
        rate ~ s(age, bs = "cs", k = k),  # k controls smoothness
        data = group_df,
        method = "REML"  # Robust smoothing parameter selection
    )
    
    return(list(fit = gam_fit, age_grid = age_grid, df = group_df))
}

convert_to_tibble <- function(smoothed) {
    # Return smoothed predictions
    tibble(
        age = smoothed[["age_grid"]],
        rate = predict(smoothed[["fit"]], data.frame(age = smoothed[["age_grid"]])),
        from = first(smoothed[["df"]]$from),
        to = first(smoothed[["df"]]$to),
        model = first(smoothed[["df"]]$model),
        transition_type = first(smoothed[["df"]]$transition_type)
    )
}

rates_smoothed <- 
    grouped_data %>% map_df(~({
    .x %>% smooth_group_spline() %>% convert_to_tibble()
    }))

# Background mortality rates
mort_pre <- read_rds(here::here("results/hmd-mortality/hp_fit_2012.rds"))
mort_post <- read_rds(here::here("results/hmd-mortality/hp_fit_2019.rds"))

aa = 22
period = "pre"

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
