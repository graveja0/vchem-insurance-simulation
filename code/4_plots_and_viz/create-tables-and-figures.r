source(here::here("code/3_estimation_and_output/run-markov.r"))
ages <- params$ages_trace
ages <- ages[-length(ages)]

source(here::here("code/helper_functions/shared-functions.r"))
source(here::here("code/0_setup/define-parameters.r"))

calibrated <- read_rds(here::here("results/calibrated-rates/calibrated-rates-with-uncertainty.rds"))
rates_smoothed <- calibrated$point_estimates
rates_smoothed %>% plot_calibration()

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

get_occupancy <- function(res,label = "") {
    res %>% map(~({
        cbind.data.frame(as.vector(.x) ,sort(rep(ages,length(params$v_tr_names)))) %>% 
            set_names(c("occupancy","age")) %>% 
            filter(age == params$age0) %>% 
            mutate(category = params$v_tr_names)
    })) %>% 
        set_names(c("Pre","Post")) %>% 
        bind_rows(.id = "period")  %>% 
        mutate(type = label) %>% 
        spread(period,occupancy) %>% 
        mutate(change = Post - Pre) %>% 
        mutate(initial = s0_post) %>% 
        mutate(overall_change = weighted.mean(change, w = initial))
}

res_ins %>% get_occupancy("Insured") 
res_employer %>% get_occupancy("Employer")
res_nongroup %>% get_occupancy("Non-Group")
res_public %>% get_occupancy("Public")


