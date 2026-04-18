rm(list = ls())
library(truncdist)
library(VGAM)
library(latex2exp)
l <- 1; u <- 2 # search region
#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02
eta <- 0.02 # true value of eta in f

eps <- 1e-3 # 1 - mass of the signal region
alpha <- 4 # candidate value for the parameter in the benchmark model q_\alpha

# signal density and CDF
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig) # pdf
Fs <- function(x) ptrunc(x, a = l, b = u, spec = 'norm', mean = mean_sig, sd = sd_sig) # CDF

# Figuring out (mu_s - d, mu_s + d) that covers an area of 1-eps:
find_d <- function(d)
{
  pl <- Fs(mean_sig-d)
  pu <- Fs(mean_sig+d)
  return(pu-pl-1+eps)
}
sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))
r <- sol$root
M_lower <- mean_sig - r # mu_s - d
M_upper <- mean_sig + r # mu_s + d

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                         rate = bkg_rate, shape = bkg_shape)
# data generating density f
f <- function(x) eta*fs(x)+(1-eta)*fb(x)

# benchmark model density:
q <- function(x){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = alpha)
}

mean1_in_g <- 0.5*M_lower + 0.5*mean_sig # location of the first Gaussian component in g 
mean2_in_g <- 0.4*M_upper + 0.6*mean_sig # location of the second Gaussian component in g
sd_in_g <- 4*sd_sig # scale parameter for the Gaussian components in g

# proposal background with dominating component on top of the benchmark model
g <- function(x, lambda) {
  phi1 <- dtrunc(x, mean = mean1_in_g, sd = sd_in_g,
                 a = l, b = u,
                 spec = 'norm')
  phi2 <-  dtrunc(x, mean = mean2_in_g, sd = sd_in_g,
                  a = l, b = u,
                  spec = 'norm')
  q <- dtrunc(x, spec = 'pareto', a = l, b = u,
              scale = l, shape = alpha)
  
  return(lambda*(phi1+phi2) + (1-2*lambda)*q)
}
# Creating the plot: Sensitivity analysis with lambda = 0.01, 0.02, 0.03
op <- par(no.readonly = TRUE)
par(mar = c(5, 5, 4, 2))
curve(g(x, 0), from = 1.1, to = 1.6, col = 'blue',
      lwd = 4, lty = 5,
      ylab = 'Density', xlab = 'x',
      cex.axis = 2, cex.lab = 2)
curve(g(x, 0.01), add = TRUE, col = 'cyan3', lwd = 4, lty = 6,
      ylab = '', xlab = 'x')
curve(g(x, 0.02), add = TRUE, col = 'orange', lwd = 4, lty = 2,
      ylab = '', xlab = 'x')
curve(g(x, 0.03), add = TRUE, col = 'purple', lwd = 4, lty = 4,
      ylab = '', xlab = 'x')
# creating the legend
legend(x = 1.26, y = 2.6, legend = TeX('$g_{\\beta}(x; \\lambda = 0.03)$'), 
       col = 'purple', lwd = 4, lty = 4,
       bty = 'n', cex = 2,
       seg.len = 2.2)
legend(x = 1.26, y = 2.4, legend = TeX('$g_{\\beta}(x; \\lambda = 0.02)$'), 
       col = 'orange', lwd = 4, lty = 2,
       bty = 'n', cex = 2,
       seg.len = 2.2)
legend(x = 1.26, y = 2.2, legend = TeX('$g_{\\beta}(x; \\lambda = 0.01)$'), 
       col = 'cyan3', lwd = 4, lty = 6,
       bty = 'n', cex = 2,
       seg.len = 2.2)
legend(x = 1.26, y = 2,  legend = TeX('$q(x) \\propto x^{-5}$'), 
       col = 'blue', lwd = 4, lty = 5,
       bty = 'n', cex = 2,
       seg.len = 2.2)
# highlighting the signal region
abline(v = c(M_lower, M_upper), col = ggplot2::alpha('grey', 0.5), lwd = 4, lty = 2)
rect(
  xleft   = M_lower,
  xright  = M_upper,
  ybottom = 0,
  ytop    = 3,
  col     = rgb(0, 1, 0, 0.15),
  border  = NA
)
par(op)
