### Load package for processing .mat file ###
library(rmatio)
library(mice)

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


### Hierarhical Model, pull data from all conditions together ###

# Full Hierarhical Model
# P(res | mu, rho) p(sigma | mu0, sigma0)  [p(mu0, sigma0) p(mu)] Flat Prior
nDataPoint <- sub1.data$nCond

# Cond Posterior
# P(mu, rho | res, mu0, sigma0)
# P(mu0, sgima0 | rho)

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

# P(mu, rho | res, mu0)
condLlhd <- function(muEst, rhoEst, resData, mu0, sigma0)
{
  jointLlhd <- dataLogllhd(resData, muEst, rhoEst) 
  + dnorm(rhoEst, mean = mu0, sd = sigma0, log = TRUE)
  return(jointLlhd)
}

# P(mu0 | rho, sigma)
# Normal distribution with var = sigma0^2/N
condSampleMu <- function(rhoEsts, sigma0)
{
  sigmaHat <- sqrt(sigma0^2 / length(rhoEsts))
  return(rnorm(1, mean = mean(rhoEsts), sd = sigmaHat))
}

# p(sigma | rho, mu0)
# Inverse Chi-square Distribution
condSampleSigma <- function(rhoEsts, mu0)
{
  sigma2Hat    <- 1 / (length(rhoEsts) - 1) * sum((rhoEsts - mu0)^2)
  sigma2Sample <- 1 / rchisq(1, df = length(rhoEsts) - 1, ncp = sigma2Hat)
  return(sqrt(sigma2Sample))
}

### Gibbs Sampler with Nested Metropolis-Hastings ###
# MH-MCMC for one run (from the cond posterior)
mhMCMC <- function(resData, mu0, sigma0, paraInit)
{
  jumpSigma2 <- 0.004
  
  # Tuning parameter 
  propCovar <- matrix(c(jumpSigma2, 0, 0, jumpSigma2), nrow = 2, ncol = 2)
  # Gaussian (symetrical) proposal distribution, no Hasitings step required
  proposal <- function(mu0, rho0) rmvnorm(1, mean = c(mu0, rho0), sigma = propCovar)
  
  # Argumented Likelihood
  sampleLlhd <- function(para)
  {
    logValue <- tryCatch({condLlhd(para[1], para[2], resData, mu0, sigma0)}, 
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
gibbsSampler <- function(dataMatrix, nSample, mu0Init, sigma0Init)
{
  mu0Sample    <- array(0, nSample)
  sigma0Sample <- array(0, nSample)
  muSample  <- matrix(0, nrow = nDataPoint, ncol = nSample)
  rhoSample <- matrix(0, nrow = nDataPoint, ncol = nSample)
  
  # Initilization 
  mu0Hat    <- mu0Init
  sigma0Hat <- sigma0Init
  muHat  <- array(0, nDataPoint)
  rhoHat <- array(0, nDataPoint)
  
  for(i in 1 : nSample)
  {
    # p(mu_j, rho_j | data, mu0)
    for(j in 1 : nDataPoint)
    {
      paraHat <- mhMCMC(dataMatrix[j, ], mu0Hat, sigma0Hat, paraInit = c(muHat[j], rhoHat[j]))
      
      muHat[j]  <- paraHat[1]
      rhoHat[j] <- paraHat[2]
    }
    
    # p(mu0 | rho, sigma0) & p(sigma0 | rho, mu0) 
    mu0Hat    <- condSampleMu(rhoHat, sigma0Hat)
    sigma0Hat <- condSampleSigma(rhoHat, mu0Hat)
    
    # Record all the values
    mu0Sample[i]    <- mu0Hat
    sigma0Sample[i] <- sigma0Hat
    
    muSample[, i]  <- muHat
    rhoSample[, i] <- rhoHat
    
    if(i %% 100 == 0)
      cat("Percent Sampled:", i / nSample, '\n')
  }
  allSamples <- 
    list(mu0Sample = mu0Sample, sigma0Sample = sigma0Sample, muSample = muSample, rhoSample = rhoSample)
  
  return(allSamples)
}

# Gibbs Sampler for Hierarhical Model
nSample <- 2e3
posteriorHierarchical <- gibbsSampler(sub1.data$dataMatrix, nSample, mu0Init = 0.2, sigma0Init = 0.2)

mu0Sample    <- posteriorHierarchical$mu0Sample
sigma0Sample <- posteriorHierarchical$sigma0Sample

# Analysis of the Posterior Distribution
par(mfrow = c(1, 2))
# Convergence & Auto-correlation Scale
plot(1 : 2e2, mu0Sample[1:2e2], type = 'l')
acf(mu0Sample)

par(mfrow = c(1, 2))
# Convergence & Auto-correlation Scale
plot(1 : 2e2, sigma0Sample[1:2e2], type = 'l')
acf(sigma0Sample)

par(mfrow = c(1, 2))
h1 <- hist(mu0Sample,    seq(-1, 1, 0.05), freq = FALSE, xlab = 'rho0',  ylab = 'density', main = '')
lines(density(mu0Sample)) # Kernel Density Estimation
title('posterior, rho0')

h2 <- hist(sigma0Sample, seq(0, 0.04 * 40, 0.04), freq = FALSE, xlab = 'sigma0',  ylab = 'density', main = '')
lines(density(sigma0Sample)) # Kernel Density Estimation
title(('posterior, sigma0'))

# Estimates of mu and rho variable
muEst  <- array(0, sub1.data$nCond)
rhoEst <- array(0, sub1.data$nCond)
muInterval  <- matrix(0, nrow = 2, ncol = sub1.data$nCond)
rhoInterval <- matrix(0, nrow = 2, ncol = sub1.data$nCond)

# Index for 95% Confidence Interval
idxLI <- nSample * 0.025
idxUI <- nSample * 0.975
for(condJ in 1:sub1.data$nCond)
{
  muSample  <- sort(posteriorHierarchical$muSample[condJ, ])
  rhoSample <- sort(posteriorHierarchical$rhoSample[condJ, ])
  
  muEst[condJ]  <- mean(muSample)
  rhoEst[condJ] <- mean(rhoSample)
  
  muInterval[, condJ]  <- c(muSample[idxLI], muSample[idxUI])
  rhoInterval[, condJ] <- c(rhoSample[idxLI], rhoSample[idxUI])
}

par(mfrow = c(1, 2))
plot(sub1.data$chPercent, muEst, 'b', xlim = c(0, 1), ylim = c(-1.5, 2), xlab = 'percent chosen', ylab = 'mu')
arrows(sub1.data$chPercent, muInterval[1, ], sub1.data$chPercent, muInterval[2, ], angle = 90, code = 3, length = 0.05)

plot(sub1.data$chPercent, rhoEst, 'b', xlim = c(0, 1), ylim = c(-0.4, 1), xlab = 'percent chosen', ylab = 'rho')
arrows(sub1.data$chPercent, rhoInterval[1, ], sub1.data$chPercent, rhoInterval[2, ], angle = 90, code = 3, length = 0.05)
