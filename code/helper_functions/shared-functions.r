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

plot_calibration <- function(rates_smoothed) {
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
}
