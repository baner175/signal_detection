rm(list = ls())
library(truncdist)
library(latex2exp)

real_l <- 1; real_u <- 35 # actual search region
l <- log(real_l); u <- log(real_u) # search region on log-scale

# signal parameters
mean_sig <- 3.5; sd_sig <- sqrt(0.01*3.5^2)

# Signal density:
fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd)*exp(x))
}
# True background density:
fb <- function(x) dtrunc(x, spec = 'exp', a = l, b = u, rate = 1.4)

# loading the background-only sample to estimate the parameter in g_\beta visa MLE:
bkg_data <- read.table('Fermi_LAT_bkg_only.txt', header = TRUE)
y <- log(bkg_data$x)
m <- length(y)

######################## GAUSSIAN TAIL PROPOSAL BACKGROUND #####################
mu_in_g <- -1; sigma_factor_in_g <- 2

# negative log likelihood based on the gaussian tail model
g_gauss_model <- function(beta){
  gi <- dtrunc(y, spec = 'norm',
               mean = mu_in_g,
               sd = sqrt(sigma_factor_in_g*beta),
               a = l, b = u)
  return(-sum(log(gi)))
}
# MLE of beta:
bth_g_gauss <- nlminb(start = 0.01,
                      objective = g_gauss_model,
                      upper = 10, lower = 0)$par

# Gaussian tail proposal background evaluated at MLE of beta
g_gauss <- function(x) dtrunc(x, spec = 'norm',
                              mean = mu_in_g,
                              sd = sqrt(sigma_factor_in_g*bth_g_gauss),
                              a = l, b = u)
g_gauss <- Vectorize(g_gauss)

######################## EXPONENTIAL PROPOSAL BACKGROUND #####################
# negative log likelihood based on the exponential model
g_exp_model <- function(beta){
  gi <- dtrunc(y, spec = 'exp',
               rate = beta,
               a = l, b = u)
  return(-sum(log(gi)))
}
# MLE of beta:
bth_g_exp <- nlminb(start = 0.01,
                    objective = g_exp_model,
                    upper = 10, lower = 0)$par
# Exponential proposal background evaluated at MLE of beta
g_exp <- function(y) dtrunc(y, spec = 'exp',
                            rate = bth_g_exp,
                            a = l, b = u)
g_exp <- Vectorize(g_exp)

########################### PLOTTING THE DENSITIES #############################
op <- par(no.readonly = TRUE)
par(mar = c(5,5,3,2),
    mgp = c(2.5,1,0))
# histogram of the background-only sample
hist(y ,probability = TRUE, breaks = 50,
     main = '', ylab = 'Density',
     xlab = 'log(x)',
     col = 'white',
     cex.lab = 2,
     cex.axis = 2,
     xlim = c(0,3.55))
curve(fb, add = TRUE, lwd = 4, lty = 1, col = 'black') # background truth
curve(g_exp, add = TRUE, lwd = 4, lty = 5, col = 'blue') # Exponential model at MLE
curve(g_gauss, add = TRUE, lwd = 4, lty = 2, col = 'brown') # Gaussian tail model at MLE 
curve(dunif(x, l, u), add = TRUE, lwd = 4, lty = 6, col = 'red') # Uniform proposal background
# generating legends
legend(x = 0.85, y = 1.25,
       legend = c(
         TeX('$f_b(x)$'),
         TeX('$g_{\\hat{\\beta}}(x) \\propto \\exp(-\\hat{\\beta}x)$'),
         TeX('$g_{\\hat{\\beta}}(x) \\propto \\exp\\left{-\\frac{(x+1)^2}{4\\hat{\\beta}}\\right}$'),
         TeX('$g(x) \\propto 1$')
       ),
       bty = 'n',
       cex = 1.7,
       lty = c(1,5,2,6),
       lwd = 4,
       col = c('black', 'blue', 'brown', 'red'),
       y.intersp = 1,
       seg.len = 3
)
par(op)