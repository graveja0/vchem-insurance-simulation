params_ <- list(
    v_tx_names = c("2012","2018"),    # Cohorts
    n_tx = 2,    
    
    v_tr_names = c("Employer","OthPrivate","Public","Uninsured"),
    v_ab_names = c("Medicare","Death"),
    n_states = 6,
    
    horizon =  17,
    
    # mort_coef12 = mortality[["2012"]]$coefficients,
    # mort_coef18 = mortality[["2018"]]$coefficients,
    # mort_fn = MortalityLaws::HP, 
    
    Delta_t = 1, # cycle duration (1 = year)
    age0 = 0  # initial age
)

params <- 
    with(params_,{
        modifyList(params_,list(
            v_names_states = c(v_tr_names, v_ab_names), # state names
            omega = horizon/Delta_t,  # Total number of cycles
            ages = (0:(horizon/Delta_t))*Delta_t + age0, 
            alpha = length(v_ab_names),
            tau = length(v_tr_names), 
            s = length(v_tr_names)*horizon/Delta_t + length(v_ab_names) #total number of states;s=τω+α
        ))
    })

params$ages_trace <- params$ages
params$ages <- params$ages[-length(params$ages)] 

