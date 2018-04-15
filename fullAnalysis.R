### Load package for processing .mat file ###
install.packages('mice')
install.packages('rmatio')

library(rmatio)
library(mice)

install.packages('rmutil')
install.packages('pracma')
install.packages('mvtnorm')

# Use package rmutil/pracma for two-dimensional numerical integration
# Use package mvtnorm for two-dimensional Gaussian Distribution
library(rmutil)
library(pracma)
library(mvtnorm)


### Experiment Data Preprocessing ###
# Point the working dir to where the data file is located
setwd("~/Desktop/")
filename <- system.file('doublePass.mat', package='rmatio')
allData  <- read.mat('./doublePass.mat')

# Extract & preprocessing data 
expCond  <- allData$cmpXsrt
subject1 <- allData$DATAJDB
nTrial   <- 1e2

# Choose one Psychometric Curve
condJ <- 3
sub1.data <- squeeze(subject1[condJ, , ])
cmpXunq   <- sort(unique(expCond[condJ, , 1]))

# Extract Data for Each Condition
chPercent  <- array(0, allData$ncmp)
dataMatrix <- matrix(0, nrow = allData$ncmp, ncol = 4)

# Four type of Choice: ++, --, -+, +-
typeMatrix <- matrix(c(1, 1, 0, 0, 0, 1, 1, 0), nrow = 4, ncol = 2, byrow = TRUE)
for(c in 1:allData$ncmp)
{
  trialData <- sub1.data[expCond[condJ, , 1] == cmpXunq[c], ]
  chPercent[c] <- sum(trialData) / length(trialData)
  
  # Count for 4 choice type
  choiceCount <- array(0, 4) 
  countFunc   <- function(c1, c2) {sum(trialData[, 1] == c1 & trialData[, 2] == c2)}
  for(type in 1:4)
  {
    label <- typeMatrix[type, ]
    choiceCount[type] <- countFunc(label[1], label[2])
  }
  # Summarize Data
  dataMatrix[c, ] <- choiceCount
}

par(mfrow = c(1,1))
plot(cmpXunq, chPercent, type = 'b', xlab = 'Disparity', ylab = 'Percent Chosen')

# Data for Subject 1, Psychometric Curve 3
sub1.data <- list(compCond = cmpXunq, chPercent = chPercent, nCond = allData$ncmp,
                  dataMatrix = dataMatrix, typeMatrix = typeMatrix)


### Likelihood Function for Multinomial Choice Data ###
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
  probPN <- (1 - (probPP + probNN)) / 2                  # p(+-) == p(-+)
  
  # Multinomial distribution density function 
  # ResData: # responses in each of the four categories (++, --, -+, +-)
  return(dmultinom(resData, prob = c(probPP, probNN, probPN, probPN), log = TRUE))
}


### ML Estimator ###
mlEstimator <- function(resData) optim(c(0, 0), function(para) 
  -logllhd(resData, para[1], para[2]), method = "L-BFGS-B", lower = c(-Inf, -1), upper = c(Inf, 1))

MLEs <- matrix(0, nrow = sub1.data$nCond, ncol = 2)
for(condJ in 1:sub1.data$nCond)
  MLEs[condJ, ] <- mlEstimator(sub1.data$dataMatrix[condJ, ])$par


### Bayesian Estimate, Full Posterior Distribution ### 
# Gird Method
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
  levels <- c(-5, -4, -3, -2, -1.5, -1, -0.5, -0.2, -0.1, -0.05) * 1e2
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

### Full Bayesian Analysis of Three Example Conditions ###
plotPosterior <- function(resData)
{
  gridSize <- 2*1e2
  
  # mu [-0.9, 0.9]; rho [-0.9, 0.9]
  muRange  <- ppoints(gridSize) * 2 - 1
  rhoRange <- ppoints(gridSize) * 1.8 - 0.9
  
  pstDist <- llhdGrid(resData, gridSize, muRange, rhoRange) 
  pstProb <- plotDnst(pstDist, gridSize, muRange, rhoRange)
  
  return(pstProb)
}

# Plot the Full Posterior for Condition 3, 5, 7 (as example)
condJ <- 3
postDnst.cond3 <- plotPosterior(sub1.data$dataMatrix[condJ, ])
condJ <- 5
postDnst.cond5 <- plotPosterior(sub1.data$dataMatrix[condJ, ])
condJ <- 7
postDnst.cond7 <- plotPosterior(sub1.data$dataMatrix[condJ, ])


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
  hist(samples[1, ], breaks = breaks, freq = FALSE, xlab = 'mu',  ylab = 'density', main = 'histograms of mu')
  hist(samples[2, ], breaks = breaks, freq = FALSE, xlab = 'rho', ylab = 'density', main = 'histograms of rho')
}


### Analysis based on MCMC 
# Number of samples & Tuning parameter
nSample    <- 2*1e4
jumpSigma2 <- 0.01

# Sample for two conditions 
condJ    <- 3
paraInit <- MLEs[condJ, ]
postSample.cond3 <- mhMCMC(sub1.data$dataMatrix[condJ, ], nSample, paraInit, jumpSigma2)
postSample.cond3 <- postSample.cond3[, 1 : nSample %% 10 == 0]

condJ    <- 7
paraInit <- MLEs[condJ, ]
postSample.cond7 <- mhMCMC(sub1.data$dataMatrix[condJ, ], nSample, paraInit, jumpSigma2)
postSample.cond7 <- postSample.cond7[, 1 : nSample %% 10 == 0]

# Plot the comparision
par(mfrow = c(2,2))
hist(postSample.cond3[1, ], seq(-1, 1, 0.04), freq = FALSE, xlab = 'mu',  ylab = 'density', main = '')
lines(postDnst.cond3$muRange, postDnst.cond3$probMu * 1e2, type = 'l')
title('mu estimate, condition 3')

hist(postSample.cond3[2, ], seq(-1, 1, 0.04), freq = FALSE, xlab = 'rho',  ylab = 'density', main = '')
lines(postDnst.cond3$rhoRange, postDnst.cond3$probRho * 1e2, type = 'l')
title('rho estimate, condition 3')

hist(postSample.cond7[1, ], seq(-1, 1, 0.04), freq = FALSE, xlab = 'mu',  ylab = 'density', main = '')
lines(postDnst.cond7$muRange, postDnst.cond7$probMu * 1e2, type = 'l')
title('mu estimate, condition 7')

hist(postSample.cond7[2, ], seq(-1, 1, 0.04), freq = FALSE, xlab = 'rho',  ylab = 'density', main = '')
lines(postDnst.cond7$rhoRange, postDnst.cond7$probRho * 1e2, type = 'l')
title('rho estimate, condition 7')


### Analysis on each condition with MCMC Sampler, compute mean and 95% posterior interval ###
nSample    <- 5*1e4
jumpSigma2 <- 0.01

# Estimates of mu and rho variable
muEst  <- array(0, sub1.data$nCond)
rhoEst <- array(0, sub1.data$nCond)
muInterval  <- matrix(0, nrow = 2, ncol = sub1.data$nCond)
rhoInterval <- matrix(0, nrow = 2, ncol = sub1.data$nCond)

# Index for 95% Confidence Interval
idxLI <- nSample / 10 * 0.025
idxUI <- nSample / 10 * 0.975

for(condJ in 1:sub1.data$nCond)
{
  # Sample with Metropolis-Hastings MCMC
  paraInit   <- MLEs[condJ, ]
  postSample <- mhMCMC(sub1.data$dataMatrix[condJ, ], nSample, paraInit, jumpSigma2)
  
  # Thinning by taking every 10 sample
  postSample <- postSample[, 1 : nSample %% 10 == 0]
  
  muSample  <- sort(postSample[1, ])
  rhoSample <- sort(postSample[2, ])
  
  muEst[condJ]  <- mean(muSample)
  rhoEst[condJ] <- mean(rhoSample)
  
  muInterval[, condJ]  <- c(muSample[idxLI], muSample[idxUI])
  rhoInterval[, condJ] <- c(rhoSample[idxLI], rhoSample[idxUI])
}

par(mfrow = c(1, 2))
plot(sub1.data$chPercent, muEst, 'b', xlim = c(0, 1), ylim = c(-1.5, 2), xlab = 'percent chosen', ylab = 'mu')
arrows(sub1.data$chPercent, muInterval[1, ], sub1.data$chPercent, muInterval[2, ], angle = 90, code = 3, length = 0.05)
title(main= "d' estimate")

plot(sub1.data$chPercent, rhoEst, 'b', xlim = c(0, 1), ylim = c(-0.4, 1), xlab = 'percent chosen', ylab = 'rho')
arrows(sub1.data$chPercent, rhoInterval[1, ], sub1.data$chPercent, rhoInterval[2, ], angle = 90, code = 3, length = 0.05)
title(main= "rho estimate")


### Hierarhical Model, pull data from all conditions together ###

# Full Hierarhical Model
# P(res | mu, sigma) p(sigma | mu0) p(mu0) p(mu)
# Hyper-parameter
nDataPoint <- sub1.data$nCond
sigmaN     <- 0.1 # Hyperparameter, control the degree of shrinkage

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
  jumpSigma2 <- 0.008
  
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
nSample <- 5e3
posteriorHierarchical <- gibbsSampler(sub1.data$dataMatrix, nSample, mu0Init = 0.2)
mu0Sample <- posteriorHierarchical$mu0Sample

par(mfrow = c(1, 2))
# Convergence & Auto-correlation Scale
plot(1 : 2e2, mu0Sample[1:2e2], type = 'l')
acf(mu0Sample)
