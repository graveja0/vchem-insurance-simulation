if (!exists("manifest")) 
    source(here::here("code/0_setup/manifest.r"))

params <- list(
    v_tx_names = c("pre","post"),    # Cohorts
    n_tx = 2,    
    
    v_tr_names = c("Employer","OthPrivate","Public","Uninsured"),
    v_ab_names = c("Medicare","Death"),
    n_states = 6,
    
    horizon =  110,

    Delta_t = 1, # cycle duration (1 = year)
    age0 = 0  # initial age
)

params <- 
    with(params,{
        modifyList(params,list(
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

####################
# Initial Occupancy
####################

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


# Background mortality rates
mort_pre <- read_rds(here::here("results/hmd-mortality/hp_fit_2012.rds"))
mort_post <- read_rds(here::here("results/hmd-mortality/hp_fit_2019.rds"))

# Transition Intensities by type and age
calibrated <- read_rds(here::here("results/calibrated-rates/calibrated-rates-with-uncertainty.rds"))
rates_smoothed <- calibrated$point_estimates

# Create a list object with transtiion intensity matrices for each age in the model.
# Each age has two items: 2012-2013 and 2019-2021 transition intenstities (pre and post)

construct_Q_medicare <- function(aa, period = "pre") {
    if (aa < 64) {
        tmp <- construct_Q(aa, period = period)
        newnames = list(c(rownames(tmp),"Medicare"), c(colnames(tmp),"Medicare"))
        tmp <- cbind(rbind(tmp,rep(0,5)),rep(0,6))
        dimnames(tmp) = newnames
        tmp <- tmp[c(params$v_tr_names,params$v_ab_names),c(params$v_tr_names,params$v_ab_names)]
    } else if (aa == 64) {
        tmp <- construct_Q(aa, period = period)
        newnames = list(c(rownames(tmp),"Medicare"), c(colnames(tmp),"Medicare"))
        tmp <- cbind(rbind(tmp,rep(0,5)),rep(0,6))
        dimnames(tmp) = newnames
        tmp <- tmp[c(params$v_tr_names,params$v_ab_names),c(params$v_tr_names,params$v_ab_names)] 
        diag(tmp) = 0
        tmp[, c("Employer")] <- tmp[, c("OthPrivate")] <- tmp[, c("Public")] <- tmp[, c("Uninsured")] <- rep(0,6)
        tmp[c("Employer"), ] <- c(0,0,0,0,1000,tmp["Employer","Death"])
        tmp[c("OthPrivate"), ] <- c(0,0,0,0,1000,tmp["OthPrivate","Death"])
        tmp[c("Public"), ] <- c(0,0,0,0,1000,tmp["Public","Death"])
        tmp[c("Uninsured"), ] <- c(0,0,0,0,1000,tmp["Uninsured","Death"])
        diag(tmp) <- -rowSums(tmp)
    } else if (aa >= 65) {
        tmp <- construct_Q(64, period = period)
        newnames = list(c(rownames(tmp),"Medicare"), c(colnames(tmp),"Medicare"))
        tmp <- cbind(rbind(tmp,rep(0,5)),rep(0,6))
        dimnames(tmp) = newnames
        tmp <- tmp[c(params$v_tr_names,params$v_ab_names),c(params$v_tr_names,params$v_ab_names)] 
        diag(tmp) = 0
        tmp[, c("Employer")] <- tmp[, c("OthPrivate")] <- tmp[, c("Public")] <- tmp[, c("Uninsured")] <- rep(0,6)
        tmp[, c("Employer")] <- tmp[, c("OthPrivate")] <- tmp[, c("Public")] <- tmp[, c("Uninsured")] <- rep(0,6)
        if (period=="pre") {
            tmp[c("Employer"), ] <- c(0,0,0,0,1000,mort_pre(aa))
            tmp[c("OthPrivate"), ] <- c(0,0,0,0,1000,mort_pre(aa))
            tmp[c("Public"), ] <- c(0,0,0,0,1000,mort_pre(aa))
            tmp[c("Uninsured"), ] <- c(0,0,0,0,1000,mort_pre(aa))
        } else if (period=="post") {
            tmp[c("Employer"), ] <- c(0,0,0,0,1000,mort_post(aa))
            tmp[c("OthPrivate"), ] <- c(0,0,0,0,1000,mort_post(aa))
            tmp[c("Public"), ] <- c(0,0,0,0,1000,mort_post(aa))
            tmp[c("Uninsured"), ] <- c(0,0,0,0,1000,mort_post(aa))
        }
        
        diag(tmp) <- -rowSums(tmp)
    }
    return(tmp)
    
}


params <- with(params,modifyList(params,list(
    m_R = ages_trace[-length(ages_trace)] %>% map( ~ ({
        list(pre = construct_Q_medicare(.x, "pre"), 
             post = construct_Q_medicare(.x, "post"))
    }))
)))

# lifetime occupancy approach requires us to transpose transition intensity matrices 
params <- with(params,modifyList(params,list(
    m_R_t = m_R %>% map(~({
        tmp <- .x
        tmp %>% map(~(t(.x)))
    }))
)))

# this transpose command is not on the matrix, but on the list object; it inverts the
# nested list
params <- with(params, modifyList(params, list(m_R_ = m_R_t %>% transpose())))
params$m_R = params$m_R_

params <- with(params,modifyList(params,list(
    m_V = m_R %>% map(~({
        R = .x
        R %>% map(~({
            m <- .x[v_tr_names,v_tr_names] 
        }))
        
    })),
    
    m_Q = m_R %>% map(~({
        R = .x 
        R %>% map(~({
            V = .x[v_tr_names,v_tr_names]
            S = .x[v_ab_names,v_tr_names]
            zero_ <- matrix(0, nrow = length(v_tr_names)+length(v_ab_names), ncol = length(v_ab_names))
            tmp <- cbind(rbind(V,S),zero_)
            dimnames(tmp) <- list(c(v_tr_names,v_ab_names),c(v_tr_names,v_ab_names))
            tmp
        }))
    }))    
)))

# clear out original m_P
params$m_P = NULL

params <- with(params,modifyList(params,list(
    m_P = m_Q %>% map(~({
        Q = .x
        Q %>% map(~(expm(.x * Delta_t)))
    }))
)))

params <- with(params,modifyList(params,list(
    m_U = m_P %>% map(~({
        P <- .x 
        P %>% map(~(.x[v_tr_names,v_tr_names]))
    })),
    m_M = m_P %>% map(~({
        P = .x
        P %>% map(~(.x[v_ab_names,v_tr_names]))
        
    }))
)))

params <- with(params,modifyList(params,list(
    D = {
        # Create diagonal age advancement matrix
        D <- matrix(0, omega, omega)
        vec <- rep(1, omega-1)
        D[row(D)-1 == col(D)] <- vec
        D[omega,omega] = 1
        D
    }
)))

vec <-  # a simple function to return the vec of an array
    function(x) {
        y <- c(x)
        return(y)
    }

vecperm <- 
    # vecperm
    # function to calculate the vec permutation matrix K of index m,n
    # let X be a m x n matrix, and X' the transpose of X
    # then K satisfies 
    # vec(X') = K*vec(X)
    function(m, n) {
        K <- matrix(0, m * n, m * n)
        a <- matrix(0, m, n)
        
        for (i in 1:m) {
            for (j in 1:n) {
                e <- a
                e[i, j] <- 1
                K <- K + kronecker(e, t(e))
            }
        }
        
        return(K)
    }


params <- with(params,modifyList(params,list(
    bbD_ = kronecker(diag(tau), D),
    bbU_ =  m_U %>% 
        map(~(bdiag(.x))),
    K = vecperm(tau, omega)
)))

params <- with(params,modifyList(params,list(
    mUtilde = bbU_ %>% map( ~ ({
        t(K) %*% bbD_ %*% K %*% .x
    }))
)))

params <- with(params,modifyList(params,list(
    mMtilde = m_M %>% map(~({
        do.call(cbind,.x) 
    }))  
)))

params <- with(params,modifyList(params,list(
    mPtilde =  map2(mUtilde, mMtilde,  ~ ({
        rbind(cbind(.x, matrix(0, tau * omega, alpha)) ,
              cbind(.y, diag(alpha)))
    }))
)))

H = with(params,matrix(1,nrow=tau, ncol=omega))

with(params,{
    V_INS <<- v_tx_names %>% map(~({
        v_ <- matrix(c(rep(1,omega),
                       rep(1,omega),
                       rep(1,omega),
                       rep(0,omega)),nrow=tau, ncol = omega,byrow = TRUE)  
    })) %>% 
        set_names(v_tx_names)
})
with(params,{
    V_UNIN <<- v_tx_names %>% map(~({
        v_ <- matrix(c(rep(0,omega),
                       rep(0,omega),
                       rep(0,omega),
                       rep(1,omega)),nrow=tau, ncol = omega,byrow = TRUE) 
    })) %>% 
        set_names(v_tx_names)
})
with(params,{
    V_PRIVATE <<- v_tx_names %>% map(~({
        v_ <- matrix(c(rep(1,omega),
                       rep(1,omega),
                       rep(0,omega),
                       rep(0,omega)),nrow=tau, ncol = omega,byrow = TRUE) 
    })) %>% 
        set_names(v_tx_names)
})
with(params,{
    V_PUB <<- v_tx_names %>% map(~({
        v_ <- matrix(c(rep(0,omega),
                       rep(0,omega),
                       rep(1,omega),
                       rep(0,omega)),nrow=tau, ncol = omega,byrow = TRUE) 
    })) %>% 
        set_names(v_tx_names)
})

with(params,{
    V_EMPLOYER <<- v_tx_names %>% map(~({
        v_ <- matrix(c(rep(1,omega),
                       rep(0,omega),
                       rep(0,omega),
                       rep(0,omega)),nrow=tau, ncol = omega,byrow = TRUE) 
    })) %>% 
        set_names(v_tx_names)
})
with(params,{
    V_NONGROUP <<- v_tx_names %>% map(~({
        v_ <- matrix(c(rep(0,omega),
                       rep(1,omega),
                       rep(0,omega),
                       rep(0,omega)),nrow=tau, ncol = omega,byrow = TRUE) 
    })) %>% 
        set_names(v_tx_names)
})

mcwr_occupancy <- function(params, H, V) {
    with(params,{
        map(v_tx_names,~({
            U = mUtilde[[.x]]
            P = mPtilde[[.x]]
            v_ = V[[.x]]
            N = solve(diag(tau*omega)-U)
            h = vec(H) %>% as.matrix()
            not_h = 1-h
            v <- vec(v_) %>% as.matrix()
            B1 <- h %*% t(v) + 0.5 * (not_h %*% t(v)) + 0.5 * (v %*% t(not_h)) # Eq. 46
            C1 = 0.5 * (rep(1,alpha) %*%  t(v)) # Eq. 48
            R1 = rbind(cbind(B1, matrix(0, tau * omega, alpha)) ,
                       cbind(C1, diag(alpha))) 
            R2 = R1 * R1
            R3 = R1 * R1 * R1
            Z = cbind(diag(tau*omega),matrix(0,nrow=tau*omega, ncol=alpha))
            e = rep(1,s)
            rho1_ <- t(N)%*% Z %*% t(P * R1) %*% e
            rho1_
        }))
    })
}
