rm(list = ls())

library(latex2exp)

# parameters to retrieve simulation results
B <- 1e5
bkg_to_phys_ratio <- 2
n_seq <- c(50, 250, 5e2, 1e3, 2e3)
beta0 <- 2 # fixed value of the parameter beta

eta <- 0.03 # eta = 0 for Fig 5 (right) and change to 0.01, 0.02 or 0.03 to generate Fig 6

# generating the plots:
transp <- 0.8
op <- par(no.readonly = TRUE)
par(mgp = c(2.5, 0.8, 0))
if(eta == 0){
  plot(x = 500, y = 0, 
       xlim = c(50, 2e3),
       ylim = c(0.02, 0.06), type = 'n',
       ylab = 'P(Type I error)',
       xlab = 'Size of the physics data',
       cex.lab = 2, cex.axis = 2,
       xaxt = 'n'
  )
  axis(1, at = n_seq,
       labels = as.character(n_seq),
       cex.axis = 1.8)
  abline(h = 0.05, col = 'grey', lwd = 4, lty = 2)
}else{
  plot(x = 500, y = 0, 
       xlim = c(50, 2e3),
       ylim = c(0.05,1), 
       yaxt = 'n',
       xaxt = 'n',
       type = 'n',
       ylab = 'Power',
       xlab = 'Size of the physics data',
       cex.lab = 2.2, cex.axis = 2.2, cex.main = 2.2,
       main = TeX(sprintf('$\\eta = %f$', eta))
  )
  ax_marks <- seq(0.05, 1, 0.1)
  axis(1, at = n_seq,
       labels = as.character(n_seq),
       cex.axis = 1.8)
  axis(2, at = ax_marks,
       labels = as.character(ax_marks),
       cex.axis = 1.8)
  abline(h = ax_marks, lty = 2, lwd = 4,
         col = 'grey')
}

# reading files containing the empiricial powers
pow_res <- c()
for(n in n_seq){
  file_name <- paste0('Numerical_example_Sec3/',
                      'WBKG__',
                      'n_phys(', n,')_',
                      'beta0(',beta0,')_',
                      'eta(',eta,')_',
                      'B(',B,')_',
                      'bkg_to_phys(',bkg_to_phys_ratio,').txt')
  pow <- read.table(file_name, header = TRUE)
  pow_res <- rbind(pow_res, pow)
}

# power curve when beta estimated:
lines(x = n_seq, y = pow_res[,1], lty = 2, lwd = 4,
      col = ggplot2::alpha('blue', transp))
points(x = n_seq, y = pow_res[,1], cex = 2,
       col = ggplot2::alpha('blue', transp), lwd = 4, pch = 16)

# power curve when beta is fixed at btea0 = 2:
lines(x = n_seq, y = pow_res[,2], lty = 4, lwd = 4, 
      col = ggplot2::alpha('brown', transp))
points(x = n_seq, y = pow_res[,2], cex = 2,
       col = ggplot2::alpha('brown', transp), lwd = 4, pch = 17)

# power curve when proposal background is uniform:
lines(x = n_seq, y = pow_res[,3], lty = 3, lwd = 4,
      col = ggplot2::alpha('red', transp))
points(x = n_seq, y = pow_res[,3], cex = 2,
       col = ggplot2::alpha('red', transp), lwd = 4, pch = 18)

# Creating legends
if(eta == 0 || eta == 0.01){
  legend('topleft',
         legend = c(
           TeX("Power-law ($\\beta$-estimated)"),
           TeX("Power-law ($\\beta = 2$)"),
           "Uniform"
         ),
         col = c("blue","brown","red"),
         lty = c(2,4,3),
         pch = 16:18,
         lwd = 4,
         cex = 2,
         seg.len = 2.2,
         bty = "n")
}else{
  legend('bottomright',
         legend = c(
           TeX("Power-law ($\\beta$-estimated)"),
           TeX("Power-law ($\\beta = 2$)"),
           "Uniform"
         ),
         col = c("blue","brown","red"),
         lty = c(2,4,3),
         pch = 16:18,
         lwd = 4,
         cex = 2,
         seg.len = 2.2,
         bty = "n")
}
par(op)