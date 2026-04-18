rm(list = ls())
library(VGAM)
library(truncdist)
library(latex2exp)

real_l <- 1; real_u <- 35  # actual search region
l <- log(real_l); u <- log(real_u) # search region on log-scale

# signal parameters:
mean_sig <- 3.5; sd_sig <- sqrt(0.01*3.5^2)

eps <- 1e-3 # 1 - mass of the signal region

# Signal density:
fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd)*exp(x))
}

# Signal CDF:
Fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(ptrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd))
}

# loading physics data to estimate the parameter of the benchmark model via MLE:
dat <- read.table('Fermi_LAT_physics.txt', header = TRUE)$x
x <- log(dat)
n <- length(x)

# negative log-likelihood using the benchmark background model q_\alpha, a shifted power-law density
q_model <- function(alpha){
  q_mass <- (1/alpha)*((l+1)^(-alpha) - (u+1)^(-alpha))
  q_i <- sapply(x, function(t){
    ((t+1)^(-alpha-1))/q_mass
  })
  return(-sum(log(q_i)))
}

alpha_hat <- nlminb(start = 0.01,
                    objective = q_model,
                    upper = Inf, lower = 0)$par

q <- function(x)
{
  q_mass <- (1/alpha_hat)*((l+1)^(-alpha_hat) - (u+1)^(-alpha_hat))
  ((x+1)^(-alpha_hat-1))/q_mass
}
# Figuring out (mu_s - d, mu_s + d) that covers an area of 1-eps:
find_d <- function(d)
{
  pl <- Fs(log(mean_sig)-d, mean = mean_sig, sd = sd_sig)
  pu <- Fs(log(mean_sig)+d, mean = mean_sig, sd = sd_sig)
  return(pu-pl-1+eps)
}
sol <- uniroot(find_d, lower = 0, upper = min(log(mean_sig) - l,u - log(mean_sig)))
r <- sol$root
M_lower <- log(mean_sig) - r # mu_s - d
M_upper <- log(mean_sig) + r # mu_s + d

# constructing the proposal background g:
# means of the Gaussian components to mix with q_\alpha
mean1_in_g <- (M_lower + log(mean_sig))/2 # location of the first Gaussian component in g
mean2_in_g <- (M_upper + log(mean_sig))/2 # location of the second Gaussian component in g
sig_fs <- sqrt(integrate(function(x) {(x^2)*fs(x)}, l, u)$value - integrate(function(x) {(x)*fs(x)}, l, u)$value^2)
sig_0 <- 3*sig_fs # SD for the Gaussian components

lambda_seq <- c(0, 0.01, 0.03, 0.05, 0.07) # lambda values to perform sensitivity analysis on
######### Plotting Densities ###################################################
# defining the proposal background g_\beta with a dominating component 
# on top of the estimated benchmark density
g <- function(y, lambda){
  q <- q(y) 
  phi1 <- dtrunc(y, spec = 'norm', a = l, b = u,
                 mean = mean1_in_g, sd = sig_0) # first Gaussian component in g
  phi2 <- dtrunc(y, spec = 'norm', a = l, b = u,
                 mean = mean2_in_g, sd = sig_0) # second Gaussian component in g
  g <- lambda*(phi1+phi2) + (1-2*lambda)*q # mixing the benchmark model with the Gaussian components
  return(g)
}

# Generating the plot for the sensitivity analysis
op <- par(no.readonly = TRUE)
par(mar = c(5,5,3,2),
    mgp = c(2.5,1,0))
# histogram of the physics data
hist(x, probability = TRUE, breaks = 50,
     main = '', ylab = 'Density',
     xlab = 'log(x)',
     col = 'white',
     cex.lab = 2,
     cex.axis = 2,
     xlim = c(0,3.55))
mycols <- c('black', 'brown', 'darkgreen', 'orange', 'purple')
palette(mycols)
my_lty <- c(3,2,6,5,4)
# plotting g_{\hat\beta} with different sizes for the dominating component
for(j in 1:length(lambda_seq))
{
  curve(g(x, lambda = lambda_seq[j]),
        l, u, add = TRUE, lwd = 4,
        col = mycols[j],
        lty = my_lty[j])
}
# highlighting the signal region
abline(v = c(M_lower, M_upper), col = 'grey', lty = 2, lwd = 4)
rect(
  xleft   = M_lower,
  xright  = M_upper,
  ybottom = 0,
  ytop    = 3,
  col     = rgb(0, 1, 0, 0.15),
  border  = NA
)
# generating legends
legend(x = 0.95, y = 1.1,
       col = mycols,
       lty = my_lty, bty = 'n', lwd = 4,
       legend=c(TeX('$g_{\\hat{\\beta}}(\\lambda = 0) \\equiv q_{\\hat{\\alpha}}$'),
                TeX(sprintf('$g_{\\hat{\\beta}}(\\lambda = %f)$', lambda_seq[-1]))),
       cex = 2,
       x.intersp = 0.5,
       seg.len = 3,
       y.intersp = 1)
par(op)