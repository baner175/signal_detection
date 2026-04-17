rm(list = ls())
library(truncdist)
## Global parameters ##

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02

l <- 1; u <- 2 # search region

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

source('power_simulation_functions.R') # loading necessary functions for simulations

B <- 1e2 # change to B = 1e5 to replicate the results in the paper
n_phys <- 250 # change n_phys to 50, 100, 250, 500, 1e3 and 2e3 to replicate the results in the paper
bkg_to_phys_ratio <- 2; eta <- 0 # set eta to 0.01, 0.02 or 0.03 for power analysis
n_bkg <- n_phys*bkg_to_phys_ratio; beta0 <- 2 # fixed value of the slope parameter of the Pareto Type I density


res_wbkg_knw <- c()
res_wbkg_est <- c()
res_wbkg_unif <- c()
start_time <- Sys.time()

# power simulation in the presence of a background with Pareto Type I proposal background, rate fixed at beta0
res_wbkg_knw <- sim_power_wbkg(eta = eta, n_phys = n_phys,
                               r = bkg_to_phys_ratio, nsims = B,
                               seed = 12345, signif.level = 0.05,
                               beta0 = beta0)

# power simulation in the presence of a background with Pareto Type I proposal background, rate parameter estimated via MLE
res_wbkg_est <- sim_power_wbkg(eta = eta, n_phys = n_phys,
                               r = bkg_to_phys_ratio, nsims = B,
                               seed = 12345, signif.level = 0.05)

# power simulation in the presence of a background with uniform proposal background
res_wbkg_unif <- sim_power_wbkg_unif(eta = eta, n_phys = n_phys,
                                     r = bkg_to_phys_ratio, nsims = B,
                                     seed = 12345, signif.level = 0.05)

end_time <- Sys.time()

# Time elapsed:
end_time - start_time 

# storing the simulation results
power_res <- data.frame(res_wbkg_est, res_wbkg_knw, res_wbkg_unif)
colnames(power_res) <- c('beta_est', paste0('beta(', beta0,')'), 'uniform_bkg')

file_name <- paste0('Numerical_example_Sec3/',
                    'WBKG__',
                    'n_phys(', n_phys,')_',
                    'beta0(', beta0,')_',
                    'eta(',eta,')_',
                    'B(',B,')_',
                    'bkg_to_phys(',bkg_to_phys_ratio,').txt')

write.table(power_res, file = file_name)

