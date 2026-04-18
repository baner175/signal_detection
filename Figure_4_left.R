rm(list = ls())

library(truncdist)
library(VGAM)
library(latex2exp)

l <- 1; u <- 2 # search region

#parameters for the signal density
mean_sig <- 1.28
sd_sig <- 0.02
# signal density
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5
# true bkg density
fb <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)

# background model for the spurious signal method
beta0 <- 4
q <- function(x, beta = beta0){ # misspecified benchmark model for background
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta)
} 
g <- function(x, lambda) lambda*fs(x) + (1-lambda)*q(x) # misspecified background with spurious signal component

# function to calculate the compensator for each value of epsilon, i.e., size of spurious signal component
delta_vs_eps <- function(eps){
  S <- function(x) fs(x)/g(x, lambda = eps) - 1
  normS <- integrate(function(x) (S(x)^2)*g(x, lambda = eps), l, u)$value |> sqrt()
  
  S_dag <- function(x) S(x)/normS
  return(integrate(function(x) S_dag(x)*fb(x), l, u)$value)
}
delta_vs_eps <- Vectorize(delta_vs_eps)
eps_seq <- c(0, 0.005, 0.01) # candidate values for epsilon
delta_vals <- sapply(eps_seq, delta_vs_eps)

# generating the plot:
op <- par(no.readonly = TRUE)
par(mgp = c(2.5, 0.8, 0))
curve(fb, from = l, to = u, 
      col = 'black', lwd = 4, lty = 1, ylab = 'Density',
      cex.lab = 2, cex.axis = 2)
my_cols <- c('blue', 'orange', 'red')
my_lty <- c(2,6,4)

# curve(gb, l, u, col = 'red')
for(i in 1:length(eps_seq)){
  curve(g(x, lambda = eps_seq[i]), col = my_cols[i],
        lwd = 4, lty = my_lty[i],
        add = TRUE)
}
legend(x = 1.1, y = 4,
       legend = c(
         TeX('$f_b;$'),
         TeX('$\\tilde{g}(\\epsilon = 0) \\equiv q;$'),
         TeX(sprintf('$\\tilde{g}(\\epsilon = %.4f);$', eps_seq[-1]))
       ),
       col = c('black', my_cols),
       lty = c(1, my_lty),
       bty = 'n',
       lwd = 4,
       seg.len = 2,
       y.intersp = 1.5, cex = 2)
legend(x = 1.6, y = 4,
       legend = TeX(sprintf('$\\tilde{\\delta} = %.5f $', c(0, round(delta_vals, 3)))),
       bty = 'n',
       y.intersp = 1.5, cex = 2)
par(op)
