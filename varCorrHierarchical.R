# Use package rmutil/pracma for two-dimensional numerical integration
# Use package mvtnorm for two-dimensional Gaussian Distribution
install.packages('rmutil')
install.packages('pracma')
install.packages('mvtnorm')

library(rmutil)
library(pracma)
library(mvtnorm)

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

expData.equ   <- expDataSimulation(0, 0.6, nTrial)
expData.thres <- expDataSimulation(1, 0.6, nTrial)
expData.ext   <- expDataSimulation(1.5, 0.6, nTrial)
expData.ill   <- expDataSimulation(2, 0.6, nTrial)

# Set up the data matrix
dataMatrix <- matrix(0, nrow = 4, ncol = 4)
dataMatrix[1, ] <- expData.equ
dataMatrix[2, ] <- expData.thres
dataMatrix[3, ] <- expData.ext
dataMatrix[4, ] <- expData.ill

### ML Estimator ###
# Mode in the likelihood space
mlEstimator <- function(resData) optim(c(0, 0), function(para) -logllhd(resData, para[1], para[2]), method = "L-BFGS-B", lower = c(-Inf, -1), upper = c(Inf, 1))

paraMode.equ   <- mlEstimator(expData.equ)$par
paraMode.thres <- mlEstimator(expData.thres)$par
paraMode.ext   <- mlEstimator(expData.ext)$par
paraMode.ill   <- mlEstimator(expData.ill)$par

### Model of the Psychophysics Data ###

# Full Model of the Data 
# P(res | mu, sigma) p(sigma | mu0) p(mu0) p(mu)
nDataPoint <- 4
sigmaN <- 0.1

# Cond Posterior
# P(mu, sigma | res, mu0)
# P(mu0 | rho)

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
dataLogllhd <- function(resData, muEst, rhoEst)
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

# P(mu, sigma | res, mu0)
condLlhd <- function(muEst, rhoEst, resData, mu0)
{
  jointLlhd <- dataLogllhd(resData, muEst, rhoEst) 
            + dnorm(rhoEst, mean = mu0, sd = sigmaN, log = TRUE)
  return(jointLlhd)
}

# P(mu0 | rho)
condSample <- function(rhoEsts)
{
  sigmaHat <- sqrt(sigmaN^2 / length(rhoEsts))
  return(rnorm(1, mean = mean(rhoEsts), sd = sigmaHat))
}

### Gibbs Sampler with Nested Metropolis-Hastings ###
# MH-MCMC for one run (from the cond posterior)
mhMCMC <- function(resData, mu0, paraInit)
{
  jumpSigma2 <- 0.004
  
  # Tuning parameter 
  propCovar <- matrix(c(jumpSigma2, 0, 0, jumpSigma2), nrow = 2, ncol = 2)
  # Gaussian (symetrical) proposal distribution, no Hasitings step required
  proposal <- function(mu0, rho0) rmvnorm(1, mean = c(mu0, rho0), sigma = propCovar)
  
  # Argumented Likelihood
  sampleLlhd <- function(para)
  {
    logValue <- tryCatch({condLlhd(para[1], para[2], resData, mu0)}, 
                         error = function(e) {return(log(0))})
    return(logValue)
  }
  
  # Initilization
  muHat  <- paraInit[1]
  rhoHat <- paraInit[2]
  
  nBurnin <- 2e2
  for(i in 1 : nBurnin)
  {
    # Sample from jumping distribution 
    paraStar <- proposal(muHat, rhoHat)
    
    # Accept/Reject
    ratio <- exp(sampleLlhd(paraStar) - sampleLlhd(c(muHat, rhoHat)))
    
    if(runif(1) <= min(ratio, 1))
    {
      muHat  <- paraStar[1]
      rhoHat <- paraStar[2]
    }
  }
  
  return(c(muHat, rhoHat))
}

# Full Gibbs Sampler
nSample <- 1e3
gibbsSampler <- function(dataMatrix, nSample, mu0Init)
{
  mu0Sample <- array(0, nSample)
  muSample  <- matrix(0, nrow = nDataPoint, ncol = nSample)
  rhoSample <- matrix(0, nrow = nDataPoint, ncol = nSample)
  
  # Initilization 
  mu0Hat <- mu0Init
  muHat  <- array(0, nDataPoint)
  rhoHat <- array(0, nDataPoint)
 
  for(i in 1 : nSample)
  {
    # p(mu_j, rho_j | data, mu0)
    for(j in 1 : nDataPoint)
    {
      paraHat <- mhMCMC(dataMatrix[j, ], mu0Hat, paraInit = c(muHat[j], rhoHat[j]))
      
      muHat[j]  <- paraHat[1]
      rhoHat[j] <- paraHat[2]
    }
    
    # p(mu0 | rho)
    mu0Hat <- condSample(rhoHat)
    
    # Record all the values
    mu0Sample[i] <- mu0Hat
    muSample[, i]  <- muHat
    rhoSample[, i] <- rhoHat
    
    if(i %% 100 == 0)
      cat("Percent Sampled:", i / nSample, '\n')
  }
  allSamples <- 
    list(mu0Sample = mu0Sample, muSample = muSample, rhoSample = rhoSample)
  
  return(allSamples)
}

# Gibbs Sampler for Hierarhical Model
postHierarchical <- gibbsSampler(dataMatrix, nSample, mu0Init = 0.2)
mu0Sample <- postHierarchical$mu0Sample

par(mfrow = c(1, 2))
# Convergence
plot(1 : 2e2, mu0Sample[1:2e2], type = 'l')
# Auto-correlation Scale
acf(mu0Sample)
