rm(list = ls())
library(truncdist)

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02

l <- 1; u <- 2 # search region

# signal density and CDF
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig) # pdf
Fs <- function(x) ptrunc(x, a = l, b = u, spec = 'norm', mean = mean_sig, sd = sd_sig) # CDF

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)


source('power_simulation_functions.R') # loading necessary functions for simulation results

B <- 1e2 # change to B = 1e5 to replicate the results in the paper
n_phys <- 250 # change n_phys to 50, 100, 250, 500, 1e3 and 2e3 to replicate the results in the paper
eta <- 0 # set eta to 0.01, 0.02 or 0.03 for power analysis
eps <- 1e-3 # eps should be such that 1-eps is the mass of the signal region
lambda_seq <- c(0.00, 0.01, 0.02, 0.03) # values of lambda that determines the size of the dominating component in g_{\beta} 

# Calculating signal region:
# Figuring out (mu_s - d, mu_s + d) that covers an area of 1-eps:
find_d <- function(d)
{
  pl <- Fs(mean_sig - d)
  pu <- Fs(mean_sig + d)
  return(pu-pl-1+eps)
}
sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))
r <- sol$root
M_lower <- mean_sig - r # mu_s - d
M_upper <- mean_sig + r # mu_s + d

# PROPOSAL BACKGROUND DENSITY PARAMETERS:
mean1_in_g <- 0.5*M_lower + 0.5*mean_sig # location of the first Gaussian component in g_{\beta}
mean2_in_g <- 0.4*M_upper + 0.6*mean_sig # location of the second Gaussian component in g_{\beta}
sd_in_g <- 4*sd_sig # scale parameter used in the Gaussian components

neg_exp_log_lik <- function(beta, lambda)
{
  -integrate(function(x){
    phi1 <- dtrunc(x, mean = mean1_in_g, sd = sd_in_g,
                      a = l, b = u,
                      spec = 'norm')
    phi2 <-  dtrunc(x, mean = mean2_in_g, sd = sd_in_g,
                       a = l, b = u,
                       spec = 'norm')
    q <- dtrunc(x, spec = 'pareto', a = l, b = u,
                 scale = l, shape = beta)
    g <- lambda*(phi1+phi2) + (1-2*lambda)*q
    
    return(log(g)*fb(x))
  }, l, u)$value
}
neg_exp_log_lik <- Vectorize(neg_exp_log_lik)

# performing simulations and calculating the corresponding values of delta
pow_seq <- c(); delta_seq <- c()
for(lam in lambda_seq){
  # power simulation in the absence of a background for different values of lambda
  pow_val <- sim_power_wobkg(eta = eta, n_phys = n_phys,
                             lambda = lam,
                             mean1_in_g = mean1_in_g,
                             mean2_in_g = mean2_in_g,
                             sd_in_g = sd_in_g,
                             nsims = B, seed = 12345,
                             signif.level = 0.05)
  pow_seq <- c(pow_seq, pow_val)
  
  # calculating minimizer of KL divergence between benchmark model Q_\alpha and F_b
  alpha_star <- nlminb(start = 0.01, 
                      objective = neg_exp_log_lik,
                      lambda = lam,
                      lower = 1e-6,
                      upper = 10)$par
  normS <- integrate(function(x){
    fs <- fs(x)
    phi1 <- dtrunc(x, mean = mean1_in_g, sd = sd_in_g,
                      a = l, b = u,
                      spec = 'norm')
    phi2 <-  dtrunc(x, mean = mean2_in_g, sd = sd_in_g,
                       a = l, b = u,
                       spec = 'norm')
    q <- dtrunc(x, spec = 'pareto', a = l, b = u,
                 scale = l, shape = beta_star)
    g <- lam*(phi1+phi2) + (1-2*lam)*q
    S <- (fs/g - 1)
    return((S^2)*g)
  }, l, u)$value |> sqrt()
  
  # delta_{\beta_star}
  delta_val <- integrate(function(x){
    fs <- fs(x)
    fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                      a = l, b = u,
                      spec = 'norm')
    fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                       a = l, b = u,
                       spec = 'norm')
    qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                 scale = l, shape = beta_star)
    gb <- lam*(fs_val1+fs_val2) + (1-2*lam)*qb
    S1 <- (fs/gb-1)/normS
    return(S1*fb_true(x))
  }, l, u)$value
  delta_seq <- c(delta_seq, delta_val)
}

# storing the simulation results:
power_res <- data.frame(lambda = lambda_seq,
                        sim_power = pow_seq,
                        delta = delta_seq)
colnames(power_res) <- c('lambda', 'power', 'delta')
file_name <- paste0('Numerical_example_Sec4/',
                    'WOBKG__',
                    'eta(',eta,')_',
                    'lambda(', paste0(lambda_seq, collapse = '_'),')_',
                    'B(',B,')_',
                    'n_phys(',n_phys,').csv')
write.csv(x = power_res, file_name,
          row.names = FALSE)
