## PBFE - TIME SERIES PREDICTION
## USING EDM
## BASED ON: "EDM FOR BEGINNERS" Chang et al. 2017

rm(list=ls())

library(rEDM)
library(dplyr)
library(ggplot2)
library(forecast)


##################################
#### running through INDIV.TS ----
##################################

# 1. stepwise ARIMA
# 2. univariate EDM
# 3. multivariate EDM

load("DATA/INDIV.TS.TRANSF.RData")

# set training data (here 60% (e.g 40) of the data is used)
midpoint = round(length(unique(INDIV.TS.TRANSF$day))*0.6)
lib_point = c(1,midpoint)
pred_point = c(midpoint+1, length(unique(INDIV.TS.TRANSF$day)))

variables = levels(INDIV.TS.TRANSF$variable)
input = INDIV.TS.TRANSF

# 1. stepwise random walk and
# 2. univariate EDM
IDs = levels(INDIV.TS.TRANSF$ID)
for (i in 1:length(IDs)){
  for (j in 1:length(variables)){
    
    # single out variable and ID and make horizontal df
    temp.df <- filter(input, ID == IDs[i]) %>%
      select(day, variable, ID, scaled) %>%
      spread(key = variable, value = scaled) %>%
      ungroup() %>%
      select(-day, -ID) 
    
    

    temp.variable <- as.data.frame(select(temp.df, variables[j]))
    temp.variable$time <- c(1:65)
    temp.variable <- temp.variable[,c(2,1)]
    
    # run random walk
    fit <- Arima(as.ts(temp.variable[1:midpoint,2]), order = c(0,1,0))
    fit2 <- ets(temp.variable[midpoint:length(max(temp.variable$time)),2], model = fit)
    onestep <- fitted(fit2)
    
    # run univar edm
    
    simp.tmp <- simplex(temp.variable, E=1:10, silent = T)
    bestE <- simp.tmp[which.min(simp.tmp$mae),"E"] # estimates optimal embedding dimension
    
    # check for non-linearity
    smap <- s_map(temp.variable, E=bestE, lib=lib_point,
                  pred=pred_point, theta=seq(0,2,0.1))
    # if best rho with theta > 0, then non-linear
    best_theta <- smap[which.max(smap$rho),"theta"][1]

    
    
    univar <- simplex(temp.variable, E=bestE, lib=lib_point,
                      pred=pred_point, stats_only = F, silent = T)
    rho <- univar[[1]]$stats$rho
    
    temp <- data.frame(rho)
    temp$type <- as.factor("univar")
    temp$E <- bestE
    temp$rmse <- univar[[1]]$stats$rmse
    temp$best.theta <- best_theta
    
    
    temp$ID <- as.factor(IDs[i])
    temp$variable <- as.factor(variables[j])
    
    if (i == 1 & j == 1){
      output.univar <- temp
    }
    else{
      output.univar <- rbind(output.univar, temp)
    }
    
  }
}