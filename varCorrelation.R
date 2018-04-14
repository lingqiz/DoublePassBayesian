# Use package rmutil/pracma for two-dimensional numerical integration
# Use package mvtnorm for two-dimensional Gaussian Distribution

install.packages('rmutil')
install.packages('pracma')
install.packages('mvtnorm')

library(rmutil)
library(pracma)
library(mvtnorm)

### Model of the Psychophysics Data ###
# P(res | mu, sigma) p(mu, sigma)

# Gaussian density wrapper
densityGauss <- function(x, y, mu, rho)
{
  # var = 1, rho 
  covar <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
  input <- matrix(0, nrow = length(x), ncol = 2)
  
  input[, 1] <- x
  input[, 2] <- y
  
  return(dmvnorm(input, mean = c(mu, mu), sigma = covar))
}

# Log likelihood function with non-informative (flat) prior on mu and sigma 
logllhd <- function(resData, muEst, rhoEst)
{
  # Bound on parameter
  bndThres = 2 * 1e-2
  if(rhoEst <= -1 + bndThres | rhoEst >= 1 - bndThres)
    return(log(0))
  
  funcWrapper <- function(x, y) return (densityGauss(x, y, muEst, rhoEst))
  
  # Probability for different quadrant
  infBnd <- 2*1e1
  probPP <- integral2(funcWrapper, xmin = 0, xmax = infBnd, ymin = 0, ymax = infBnd, 
                      reltol = 1e-2, maxlist = 2*1e4)$Q  # p(++)
  probNN <- integral2(funcWrapper, xmin = -infBnd, xmax = 0, ymin = -infBnd, ymax = 0, 
                      reltol = 1e-2, maxlist = 2*1e4)$Q  # p(--)  
  probPN <- (1 - (probPP + probNN)) / 2  # p(+-) == p(-+)
  
  # Multinomial distribution density function 
  # ResData: # responses in each of the four categories (++, --, -+, +-)
  return(dmultinom(resData, prob = c(probPP, probNN, probPN, probPN), log = TRUE))
}


### Exp Data Simulation ###
# Simulation of response "data" with true parameter mu and rho
expDataSimulation <- function(mu, rho, nTrial)
{
  covar <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
  trialData <- array(0, 4) # ++, --, -+, +-
  
  # Decision variable
  dv <- rmvnorm(nTrial, mean = c(mu, mu), sigma = covar)
  
  # Response "data"
  trialData[1] = sum((dv[, 1] > 0) & (dv[, 2] > 0))
  trialData[2] = sum((dv[, 1] < 0) & (dv[, 2] < 0))
  trialData[3] = sum((dv[, 1] < 0) & (dv[, 2] > 0))
  trialData[4] = sum((dv[, 1] > 0) & (dv[, 2] < 0))
  
  return(trialData)
}

nTrial <- 1e3
expData.rnd   <- expDataSimulation(0, 0, nTrial)
expData.equ   <- expDataSimulation(0, 0.6, nTrial)
expData.thres <- expDataSimulation(1, 0.6, nTrial)
expData.ext   <- expDataSimulation(1.5, 0.6, nTrial)
expData.ill   <- expDataSimulation(2, 0.6, nTrial)

### ML Estimator ###
# Mode in the likelihood space
mlEstimator <- function(resData) optim(c(0, 0), function(para) -logllhd(resData, para[1], para[2]), method = "L-BFGS-B", lower = c(-Inf, -1), upper = c(Inf, 1))

paraMode.rnd   <- mlEstimator(expData.rnd)$par   
paraMode.equ   <- mlEstimator(expData.equ)$par
paraMode.thres <- mlEstimator(expData.thres)$par
paraMode.ext   <- mlEstimator(expData.ext)$par
paraMode.ill   <- mlEstimator(expData.ill)$par

### Bayesian Estimate, Full Posterior Distribution ### 
### Grid Sampler ###
llhdGrid <- function(resData, gridSize, muRange, rhoRange)
{
  logValue <- matrix(0, nrow = gridSize, ncol = gridSize)
  for(i in 1 : gridSize)
  {
    for(j in 1 : gridSize)
    {
      logll <- tryCatch(
        {logllhd(resData, muRange[i], rhoRange[j])},
        error = function(e)
        {
          print(c(i, j))
          if(j == 1)
            return(logValue[i - 1, 1])
          else
            return(logValue[i, j - 1])
        })
      logValue[i, j] <- logll
    }
  }
  return(logValue)
}

plotDnst <- function(llhdGrid, gridSize, muRange, rhoRange)
{
  par(mfrow = c(1,3))
  levels <- c(-15, -10, -6, -5, -4, -3, -2, -1.5, -1, -0.5, -0.2) * 1e2
  contour(muRange, rhoRange, llhdGrid, xlab = "mu", ylab = "rho", levels = levels, pch = 1)
  
  # Marginal Distribution 
  probMu  <- array(0, gridSize)
  probRho <- array(0, gridSize)
  
  for(i in 1 : gridSize)
    for(j in 1 : gridSize)
      probMu[i] <- probMu[i] + exp(llhdGrid[i, j])
  
  for(j in 1 : gridSize)
    for(i in 1 : gridSize)
      probRho[j] <- probRho[j] + exp(llhdGrid[i, j])
  
  probMu  <- probMu / sum(probMu)
  probRho <- probRho / sum(probRho)
  
  plot(muRange, probMu, type = 'l', xlab = "mu", ylab = "density", pch = 2)
  plot(rhoRange, probRho, type = 'l', xlab = "rho", ylab = "density", pch = 3)
  
  probObj <- data.frame("muRange" = muRange, "probMu" = probMu,
                        "rhoRange" = rhoRange, "probRho" = probRho)
  return(probObj)
}

### Condition 1 Rnd
gridSize <- 2*1e2

# mu [-0.9, 0.9]; rho [-0.9, 0.9]
muRange  <- ppoints(gridSize) * 1.8 - 0.9
rhoRange <- ppoints(gridSize) * 1.8 - 0.9

pstDist.rnd <- llhdGrid(expData.rnd, gridSize, muRange, rhoRange) 
prob.rnd <- plotDnst(pstDist.rnd, gridSize, muRange, rhoRange)

### Condition 2 Equ 
gridSize <- 2*1e2

# mu [-0.9, 0.9]; rho [-0.9, 0.9]
muRange  <- ppoints(gridSize) * 1.8 - 0.9
rhoRange <- ppoints(gridSize) * 1.8 - 0.9

pstDist.equ <- llhdGrid(expData.equ, gridSize, muRange, rhoRange) 
prob.equ <- plotDnst(pstDist.equ, gridSize, muRange, rhoRange)

### Condition 3 Thres 
gridSize <- 2*1e2

# mu [0, 2]; rho [-0.9, 0.9]
muRange  <- ppoints(gridSize) * 2
rhoRange <- ppoints(gridSize) * 1.8 - 0.9

pstDist.thres <- llhdGrid(expData.thres, gridSize, muRange, rhoRange) 
prob.thres <- plotDnst(pstDist.thres, gridSize, muRange, rhoRange)

### Condition 4 Extrem
gridSize <- 2*1e2

# mu [0.2, 2.2]; rho [-0.9, 0.9]
muRange  <- ppoints(gridSize) * 2 + 0.2
rhoRange <- ppoints(gridSize) * 1.8 - 0.9

pstDist.ext <- llhdGrid(expData.ext, gridSize, muRange, rhoRange) 
prob.ext <- plotDnst(pstDist.ext, gridSize, muRange, rhoRange)

### Condition 5 ill-Defined
gridSize <- 2*1e2

# mu [0.5, 2.5]; rho [-0.9, 0.9]
muRange  <- ppoints(gridSize) * 2 + 0.5
rhoRange <- ppoints(gridSize) * 1.8 - 0.9

pstDist.ill <- llhdGrid(expData.ill, gridSize, muRange, rhoRange)
prob.ill <- plotDnst(pstDist.ill, gridSize, muRange, rhoRange)

### Metropolis-Hastings Sampler ###
# MH MCMC Implementation
mhMCMC <- function(resData, nSample, paraInit, jumpSigma2)
{
  postSample <- matrix(0, nrow = 2, ncol = nSample)
  
  # Tuning parameter 
  propCovar <- matrix(c(jumpSigma2, 0, 0, jumpSigma2), nrow = 2, ncol = 2)
  # Gaussian (symetrical) proposal distribution, no Hasitings step required
  proposal <- function(mu0, rho0) rmvnorm(1, mean = c(mu0, rho0), sigma = propCovar)
  
  # Argumented Likelihood
  sampleLlhd <- function(para)
  {
    logValue <- tryCatch({logllhd(resData, para[1], para[2])}, 
                         error = function(e) {return(log(0))})
    return(logValue)
  }
  
  # Initilization
  muHat  <- paraInit[1]
  rhoHat <- paraInit[2]
  
  nAccept <- 0
  for(i in 1 : nSample)
  {
    # Sample from jumping distribution 
    paraStar <- proposal(muHat, rhoHat)
    
    # Accept/Reject
    ratio <- exp(sampleLlhd(paraStar) - sampleLlhd(c(muHat, rhoHat)))
    
    if(runif(1) <= min(ratio, 1))
    {
      nAccept <- nAccept + 1
      muHat  <- paraStar[1]
      rhoHat <- paraStar[2]
    }
      
    # muSample
    postSample[1, i] <- muHat
    # rhoSample
    postSample[2, i] <- rhoHat
  }
  
  cat('MCMC Accept Rate:', nAccept / nSample, '\n')
  return(postSample)
}

convCheck <- function(nPlot, samples)
{
  par(mfrow = c(2,2))
  
  # Convergence Check
  plot(1 : nPlot, samples[1, 1 : nPlot], type = 'l')
  plot(1 : nPlot, samples[2, 1 : nPlot], type = 'l')
  
  # Auto-correlation Scale
  acf(samples[1, 1 : nPlot])
  acf(samples[2, 1 : nPlot])
}

histPlot <- function(samples, breaks)
{
  par(mfrow = c(1, 2))
  hist(samples[1, ], breaks = breaks, freq = FALSE, xlab = 'mu',  ylab = 'density', main = 'Histograms of mu')
  hist(samples[2, ], breaks = breaks, freq = FALSE, xlab = 'rho', ylab = 'density', main = 'Histograms of rho')
}

# Number of samples 
nSample <- 2*1e4
nPlot <- 1e3

### Condition 1, Rnd
paraInit <- paraMode.rnd
postSample.rnd <- mhMCMC(expData.rnd, nSample, paraInit, jumpSigma2 = 0.004)

# Convergence Check
convCheck(nPlot, postSample.rnd)
# Thinning by taking every 10 sample
postSample.rnd <- postSample.rnd[, 1 : nSample %% 10 == 0]

# Histogram of the distribution
histPlot(postSample.rnd, breaks = 30)

### Condition 2, Equ
paraInit <- c(0, 0)
postSample.equ <- mhMCMC(expData.equ, nSample, paraInit, jumpSigma2 = 0.004)

# Convergence Check 
convCheck(nPlot, postSample.equ)
# Thinning by taking every 20 sample
postSample.equ <- postSample.equ[, 1 : nSample %% 20 == 0]

# Histogram of the distribution
histPlot(postSample.equ, breaks = 25)

### Condition 4, Ext
paraInit <- c(0, 0)
postSample.ext <- mhMCMC(expData.ext, nSample, paraInit, jumpSigma2 = 0.004)

# Convergence Check 
convCheck(nPlot, postSample.ext)
# Thinning & Discard Burning
postSample.ext <- postSample.ext[, 1 : nSample %% 20 == 0]
postSample.ext <- postSample.ext[, 10 : length(postSample.ext[1, ])]

# Histogram of the distribution
histPlot(postSample.ext, breaks = 30)

