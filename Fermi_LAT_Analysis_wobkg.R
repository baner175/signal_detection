rm(list = ls())

################################################################################
################################################################################
##                                                                            ##
## Code needed to reproduce the analysis of the simulated data from           ##
## Fermi-LAT without a background-only sample                                 ##
##                                                                            ##
################################################################################
################################################################################

####################### Loading the required packages ##########################
library(VGAM)
library(truncdist)
library(knitr)
################################################################################

# Loading the data and transforming into log-scale
dat <- read.table('Fermi_LAT_physics.txt', header = TRUE)
x <- dat$x
y <- log(x) # changing the data to log-scale
n <- length(y)

real_l <- 1; real_u <- 35 # Search region
l <- log(real_l); u <- log(real_u) # Search region on log-scale

# parameters for the signal density
mean_sig <- 3.5; sd_sig <- sqrt(0.01*3.5^2)
eps <- 1e-3 # mass outside the signal region

# SIGNAL DENSITY:
fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd)*exp(x))
}

# SIGNAL CDF:
Fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(ptrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd))
}

# Finding the signal region around log(mean_sig) with mass 1-eps:
find_d <- function(d)
{
  pl <- Fs(log(mean_sig)-d, mean = mean_sig, sd = sd_sig)
  pu <- Fs(log(mean_sig)+d, mean = mean_sig, sd = sd_sig)
  return(pu-pl-1+eps)
}

sol <- uniroot(find_d, lower = 0, upper = min(log(mean_sig) - l,u - log(mean_sig)))

r <- sol$root

M_lower <- log(mean_sig) - r # lower bound of the signal region
M_upper <- log(mean_sig) + r # upper bound of the signal region


# negative log-likelihood using the nominal background model
q_model <- function(alpha){
  q_mass <- (1/alpha)*((l+1)^(-alpha) - (u+1)^(-alpha))
  q_i <- sapply(y, function(t){
    ((t+1)^(-alpha-1))/q_mass
  })
  return(-sum(log(q_i)))
}
alpha_hat <- nlminb(start = 0.01,
                   objective = q_model,
                   upper = 10, lower = 0)$par # MLE of alpha

# defining q_alpha at MLE:
q <- function(x) {
  q_mass <- (1/alpha_hat)*((l+1)^(-alpha_hat) - (u+1)^(-alpha_hat))
  return(((x+1)^(-alpha_hat-1))/q_mass)
}

# constructing the proposal background g:
# means of the Gaussian components to mix with q_\alpha
mean1_in_g <- (M_lower + log(mean_sig))/2 
mean2_in_g <- (M_upper + log(mean_sig))/2
sig_fs <- sqrt(integrate(function(x) {(x^2)*fs(x)}, l, u)$value - integrate(function(x) {(x)*fs(x)}, l, u)$value^2)
sig_0 <- 3*sig_fs # SD for the gaussian components

# Defining proposal background g:
g <- function(x, lambda){
  phi_1 <- dtrunc(x, spec = 'norm', a = l, b = u,
                    mean = mean1_in_g, sd = sig_0)
  phi_2 <- dtrunc(x, spec = 'norm', a = l, b = u,
                    mean = mean2_in_g, sd = sig_0)
  return(
    lambda*(phi_1 + phi_2) + (1-2*lambda)*q(x)
  )
}

d_log_h <- function(t) -log(t+1)
E_q_d_log_h <- integrate(function(t) d_log_h(t)*q(t), l, u)$value
d_log_q <- function(t) d_log_h(t) - E_q_d_log_h
d2_log_q_int1 <- integrate(function(t) {
  (d_log_h(t)^2)*q(t)
}, l, u)$value
d2_log_q <- -d2_log_q_int1 + E_q_d_log_h^2
signal_search <- function(lambda){
  norm_S <- integrate(function(x){
    return(((fs(x)/g(x, lambda = lambda) - 1)^2)*g(x, lambda = lambda))
  }, l, u)$value |> sqrt()
  
  S0_vec <- sapply(y, function(x){
    S_val <- (fs(x)/g(x, lambda = lambda) - 1)
    return(S_val/(norm_S^2))
  })
  theta0_hat <- mean(S0_vec)
  d_log_q_vec <- sapply(y, d_log_q)
  d_normS_sq <- -(1-2*lambda)*integrate(function(y){
    q <- q(y)
    fs <- fs(y)
    g <- g(y, lambda = lambda)
    d_log_q <- d_log_q(y)
    return(((fs/g)^2)*q*d_log_q)
  }, l, u)$value
  
  d_S0_vec <- sapply(y, function(y){
    fs <- fs(y)
    q <- q(y)
    g <- g(y, lambda = lambda)
    d_log_q <- d_log_q(y)
    return(-((norm_S^2)*(fs/(g^2))*(1-2*lambda)*q*d_log_q + (fs/g-1)*d_normS_sq)/(norm_S^4))
  })
  V_hat <- mean(d_log_q_vec^2)
  J_hat <- -d2_log_q
  d_theta0 <- mean(d_S0_vec)
  var_S0_F_hat <- mean(S0_vec^2) - theta0_hat^2
  sig_theta0_hat <- sqrt(
    var_S0_F_hat + (V_hat/(J_hat^2))*(d_theta0^2) + 
      (2/J_hat)*d_theta0*
      mean(S0_vec*d_log_q_vec)
  )
  theta0_stat <- sqrt(n)*(theta0_hat-0)/sig_theta0_hat
  p_val <- pnorm(theta0_stat, lower.tail = FALSE)
  return(c(theta0_hat,
           p_val))
}

# lambda_max <- 0.05 # change lambda_max here
lambda_seq <- c(0.03, 0.05, 0.07)

res_sig_search <- sapply(lambda_seq, signal_search)
res_sig_search <- cbind(lambda_seq, t(res_sig_search))
res_sig_search <- data.frame(res_sig_search)
colnames(res_sig_search) <- c("lambda", "theta0beta_hat", "p-value")
caption <- paste0("Signal Search Results")
kable(res_sig_search, format = "simple", digits = 10,
      booktabs = TRUE, escape = FALSE,
      caption = caption)
