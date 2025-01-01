
if (!exists("manifest")) 
    source(here::here("code/0_setup/manifest.r"))
if (!exists("params"))
    source(here::here("code/0_setup/define-parameters.r"))

if (!file.exists(here("_data/raw/hmd/life-table-usa.rds"))) {
    source("~/auth-aws.r")
    tmp <- demography::hmd.mx("USA",username = hmd_u, password = hmd_p, paste0("USA"))
    write_rds(tmp,file = here("_data/raw/hmd/life-table-usa.rds"))
} else {
    tmp <- read_rds(here("_data/raw/hmd/life-table-usa.rds"))
}


# 2012 Mortality

lt =  demography::lifetable(tmp,series = "male", years = 2012)   # life table for background mortality

age0 = params$age0
time_step = 1/4
horizon = 111
n_cycles = horizon / time_step
ages_in_model = (0:(horizon/time_step))*time_step + age0

ages     <- lt$age
deaths   <- lt$lx * lt$qx
exposure <- lt$lx
# restrict to only those in specified age ragne
deaths <- deaths[ages>=age0]
exposure <- exposure[ages>=age0]
ages <- ages[ages>=age0]

mort_fit <- MortalityLaws::MortalityLaw(
    x  = ages,
    Dx  = deaths,   # vector with death counts
    Ex  = exposure, # vector containing exposures
    law = "HP2",
    opt.method = "LF2")
plot(mort_fit)
# Extract Death Rates by Age

r_ <- MortalityLaws::HP2(ages,mort_fit$coefficients)$hx
yrs = c(1,4,diff(ages[-c(1)]))
r_ <- r_ / yrs
afn_ <- Hmisc::approxExtrap(ages,r_,xout = ages_in_model)
afn_ <- approxfun(afn_$x,afn_$y)

# check
plot(age0:111,afn_(age0:111))
points(ages_in_model,afn_(ages_in_model),col = "red")

r <- tibble(index = 0:n_cycles,
            age = ages_in_model,
            r_die = afn_(ages_in_model)) %>%
    mutate(yrs = c(time_step,diff(age))) %>%
    mutate(r_die = r_die * yrs)

# Mini Markov Alive-Dead Model to Validate Life Expectancy
s_ <- c(1,0)
tr <- 0:n_cycles %>% map(~({
    px <- 1-exp(-r[.x+1,]$r_die)
    P = matrix(c(1-px,0,px,1),nrow=2,ncol=2, dimnames = list(c("A","D"),c("A","D")))
    s_ <<- s_ %*% P
    data.frame(s_)
})) %>%
    bind_rows() %>%
    as.matrix()
# Everyone starts alive
tr <- rbind(c(1,0),tr)

# Half-Cycle Correction
hc <- rep(1,nrow(tr))
hc[1] = 0.5
hc[length(hc)] = .5
le <- sum((tr %*% c(1,0)*time_step)*hc)

cat(paste0("Life Expectanccy (Modeled):",round(le+age0-1,2),"\n"))
cat(paste0("Life Expectancy (Life Table):", round(lt$ex[1],2)))

afn_ %>% write_rds(here::here("results/hmd-mortality/hp_fit_2012.rds"))

lt =  demography::lifetable(tmp,series = "male", years = 2019)   # life table for background mortality

age0 = params$age0
time_step = 1/4
horizon = 111
n_cycles = horizon / time_step
ages_in_model = (0:(horizon/time_step))*time_step + age0

ages     <- lt$age
deaths   <- lt$lx * lt$qx
exposure <- lt$lx
# restrict to only those in specified age ragne
deaths <- deaths[ages>=age0]
exposure <- exposure[ages>=age0]
ages <- ages[ages>=age0]

mort_fit <- MortalityLaws::MortalityLaw(
    x  = ages,
    Dx  = deaths,   # vector with death counts
    Ex  = exposure, # vector containing exposures
    law = "HP2",
    opt.method = "LF2")
plot(mort_fit)
# Extract Death Rates by Age

r_ <- MortalityLaws::HP2(ages,mort_fit$coefficients)$hx
yrs = c(1,4,diff(ages[-c(1)]))
r_ <- r_ / yrs
afn_ <- Hmisc::approxExtrap(ages,r_,xout = ages_in_model)
afn_ <- approxfun(afn_$x,afn_$y)

# check
plot(age0:111,afn_(age0:111))
points(ages_in_model,afn_(ages_in_model),col = "red")

r <- tibble(index = 0:n_cycles,
            age = ages_in_model,
            r_die = afn_(ages_in_model)) %>%
    mutate(yrs = c(time_step,diff(age))) %>%
    mutate(r_die = r_die * yrs)

# Mini Markov Alive-Dead Model to Validate Life Expectancy
s_ <- c(1,0)
tr <- 0:n_cycles %>% map(~({
    px <- 1-exp(-r[.x+1,]$r_die)
    P = matrix(c(1-px,0,px,1),nrow=2,ncol=2, dimnames = list(c("A","D"),c("A","D")))
    s_ <<- s_ %*% P
    data.frame(s_)
})) %>%
    bind_rows() %>%
    as.matrix()
# Everyone starts alive
tr <- rbind(c(1,0),tr)

# Half-Cycle Correction
hc <- rep(1,nrow(tr))
hc[1] = 0.5
hc[length(hc)] = .5
le <- sum((tr %*% c(1,0)*time_step)*hc)

cat(paste0("Life Expectanccy (Modeled):",round(le+age0-1,2),"\n"))
cat(paste0("Life Expectancy (Life Table):", round(lt$ex[1],2)))

afn_ %>% write_rds(here::here("results/hmd-mortality/hp_fit_2019.rds"))
