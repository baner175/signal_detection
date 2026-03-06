rm(list = ls())

library(VGAM)
library(truncdist)
library(latex2exp)
library(knitr)
library(kableExtra)
library(ggplot2)

real_l <- 1; real_u <- 35
l <- log(real_l); u <- log(real_u)

mean_sig <- 3.5; sd_sig <- sqrt(0.01*3.5^2)
eps <- 1e-3
mu_in_qb <- -1; sigma_factor_in_qb <- 2

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

fb <- function(x) dtrunc(x, spec = 'exp', a = l, b = u, rate = 1.4)

phys_data <- read.table('Fermi_LAT_physics.txt', header = TRUE)
bkg_data <- read.table('Fermi_LAT_bkg_only.txt', header = TRUE)
y <- log(bkg_data$x)
x <- log(phys_data$x)
n <- length(x)
m <- length(y)


g_gauss_model <- function(beta){
  gi <- dtrunc(y, spec = 'norm',
               mean = mu_in_qb,
               sd = sqrt(sigma_factor_in_qb*beta),
               a = l, b = u)
  return(-sum(log(gi)))
}

bth_g_gauss <- nlminb(start = 0.01,
                      objective = g_gauss_model,
                      upper = 10, lower = 0)$par

g_gauss <- function(x) dtrunc(x, spec = 'norm',
                              mean = mu_in_qb,
                              sd = sqrt(sigma_factor_in_qb*bth_g_gauss),
                              a = l, b = u)
g_gauss <- Vectorize(g_gauss)

g_exp_model <- function(beta){
  gi <- dtrunc(y, spec = 'exp',
               rate = beta,
               a = l, b = u)
  return(-sum(log(gi)))
}

bth_g_exp <- nlminb(start = 0.01,
                    objective = g_exp_model,
                    upper = 10, lower = 0)$par
g_exp <- function(y) dtrunc(y, spec = 'exp',
                            rate = bth_g_exp,
                            a = l, b = u)
g_exp <- Vectorize(g_exp)

par(mar = c(5,5,3,2),
    mgp = c(2.5,1,0))
hist(y ,probability = TRUE, breaks = 50,
     main = '', ylab = 'Density',
     xlab = 'log(x)',
     col = 'white',
     cex.lab = 2,
     cex.axis = 2,
     xlim = c(0,3.55))

curve(fb, add = TRUE, lwd = 4, lty = 1, col = 'black')
curve(g_exp, add = TRUE, lwd = 4, lty = 5, col = 'blue')
curve(g_gauss, add = TRUE, lwd = 4, lty = 2, col = 'brown')
curve(dunif(x, l, u), add = TRUE, lwd = 4, lty = 6, col = 'red')

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
