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


