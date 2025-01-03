source(here::here("code/0_setup/manifest.r"))
source(here::here("code/3_estimation_and_output/setup-markov.r"))

res_ins <- params %>% mcwr_occupancy(H = H, V = V_INS)
res_public <- params %>% mcwr_occupancy(H = H, V = V_PUB)
res_employer <- params %>% mcwr_occupancy(H = H, V = V_EMPLOYER )
res_private <- params %>% mcwr_occupancy(H = H, V = V_PRIVATE )
res_nongroup <- params %>% mcwr_occupancy(H = H, V = V_NONGROUP)
res_uninsured <- params %>% mcwr_occupancy(H = H, V = V_UNIN)

