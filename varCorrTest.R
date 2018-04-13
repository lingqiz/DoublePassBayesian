# Use package rmutil for two-dimensional numerical integration

install.packages('rmutil')
install.packages('pracma')
install.packages('mvtnorm')

library(rmutil)
library(pracma)
library(mvtnorm)

densityGauss <- function(x, y, mu, rho)
{
  covar <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
  input <- matrix(0, nrow = length(x), ncol = 2)
  
  input[, 1] <- x
  input[, 2] <- y
  
  return(dmvnorm(input, mean = c(mu, mu), sigma = covar))
}

funcWrapper <- function(x, y) return (densityGauss(x, y, -0.891, 0.007))


# Probability for different quadrant
# Test different numerical integration
infBnd <- 1e2
p1 <- int2(funcWrapper, a = c(0, 0), b = c(Inf, Inf), eps = 1e-4)    # p(++)
p1 <- integral2(funcWrapper, xmin = 0, xmax = infBnd, ymin = 0, ymax = infBnd, maxlist = 1e4)$Q

p2 <- int2(funcWrapper, a = c(-Inf, -Inf), b = c(0, 0), eps = 1e-4)  # p(--)
p2 <- integral2(funcWrapper, xmin = -infBnd, xmax = 0, ymin = -infBnd, ymax = 0, maxlist = 1e4)$Q  # p(--)

# p(+-) == p(-+)
p3 <- int2(funcWrapper, a = c(-Inf, 0), b = c(0, Inf), eps = 1e-4) 
p3 <- integral2(funcWrapper, xmin = -infBnd, xmax = 0, ymin = 0, ymax = infBnd)$Q

p4 <- int2(funcWrapper, a = c(0, -Inf), b = c(Inf, 0), eps = 1e-4)
p4 <- integral2(funcWrapper, xmin = 0, xmax = infBnd, ymin = -infBnd, ymax = 0)$Q

# x: count for different "category"
dmultinom(c(1, 1, 1, 1), c(p1, p2, p3, p4), log = FALSE)

