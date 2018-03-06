## PBFE - CONVERGENT CROSS MAPPING OF SPECIES INTERACTIONS
## David Inauen - 01.03.2018
## Mostly based on Chang 2016 - EDM for beginners 
## and supp. material


rm(list=ls())

library(rEDM)
library(Kendall)
library(dplyr)
library(tidyr)
library(rlist)

# exploration


load("DATA/INDIV.TS.TRANSF.RData")

IDs <- levels(INDIV.TS.TRANSF$ID)

start.time <- Sys.time()
for (ids in 1:length(IDs)){

  
  # filtering out one ID and make horizontal dfs
  temp <- filter(INDIV.TS.TRANSF, ID == IDs[ids]) %>%
    ungroup() %>% 
    select(day, variable, scaled) %>%
    spread(key = variable, value = scaled) %>%
    ungroup() %>%
    select(-day)
  
  # design sequence of library size
  maxlib <- length(temp$temperature)
  libs <- c(round(seq(0.1, 1, 0.1)*maxlib))
  perc <- c("25%","50%","75%","100%")
  
  
  ## CCM ANALYSIS testing N1 cross-mapping N2 (test N2 as causation for N1)
  
  # go through all species combinations
  # generate a df where all variable are crossmapped against each other 
  # loop always checks that no redundant calculations are performed
  
  variables.done <- list(0)
  for (n1 in 1:ncol(temp)){
    for (n2 in 1:ncol(temp)){
      
      N1 <- names(temp)[n1]
      N2 <- names(temp)[n2]
      
      if (n1 != n2){
        
        # ensures that no reduntant calculations are performed
        if (N2 %in% variables.done){
          print("We've done that!")
        }
        else{
          print(paste0("CURRENT ID", IDs[ids]))
          print(paste0("RUN:", n1))
          print(paste0("Variable N1:", N1))
          print(paste0("Variable N2:", N2))
          
          
          # determine embedding dimension
          E.test.n1 = NULL
          n12q=NULL
          for(E.t in 2:8){
            cmxy.t <- ccm(temp, E = E.t, lib_column = N1,
                          target_column = N2, lib_sizes = maxlib,
                          num_samples = 1,
                          tp=-1, random_libs = F, silent = T)
            E.test.n1 = rbind(E.test.n1,cmxy.t)}
          E_n1 <- E.test.n1$E[which.max(E.test.n1$rho)[1]]
          
          # CCM analysis with varying library size (L)
          n1_xmap_n2 <- ccm(temp, E=E_n1,lib_column= N1,
                            target_column= N2,
                            lib_sizes=libs, num_samples=200, replace=T,
                            RNGseed = 2301, silent = T)
          n1_xmap_n2 <- na.omit(n1_xmap_n2)
          n12q= as.matrix(aggregate(n1_xmap_n2[,c('rho')],
                                    by = list(as.factor(n1_xmap_n2$lib_size)),
                                    quantile)[,'x'])
          # ensures that no inf values sneak in (happens only for didinium, diss. oxygen interaction)
          n12q <- n12q[is.finite(rowSums(n12q)),]
          n1v2 <- apply(n12q[,2:5],2,MannKendall)
          
          ##################################s######################
          ## CCM ANALYSIS testing N2 cross-mapping N1 (test N1 as causation for N2)
          E.test.n2 = NULL
          n21q=NULL
          for(E.t in 2:8){
            cmxy.t <- ccm(temp, E = E.t, lib_column = N2,
                          target_column = N1, lib_sizes = maxlib,
                          num_samples = 1,
                          tp=-1, random_libs = F, silent = T)
            E.test.n2 = rbind(E.test.n2,cmxy.t)}
          E_n2 <- E.test.n2$E[which.max(E.test.n2$rho)[1]]
          
          # CCM
          n2_xmap_n1 <- ccm(temp, E=E_n2,lib_column=N2,
                            target_column=N1, 
                            lib_sizes=libs, num_samples=200, replace=T,
                            RNGseed=2301, silent = T)
          n2_xmap_n1 <- na.omit(n2_xmap_n1)
          n21q=as.matrix(aggregate(n2_xmap_n1[,c('rho')],
                                   by = list(as.factor(n2_xmap_n1$lib_size)),
                                   quantile)[,'x'])
          n21q <- n21q[is.finite(rowSums(n21q)),]
          n2v1 <- apply(n21q[,2:5],2,MannKendall)
          
          for (i in 1:4){
            pval <- as.numeric((n1v2[[i]])$sl)
            percent <- as.factor(perc[i])
            temp_out <- data.frame(pval, percent)
            temp_out$N1 <- as.factor(N1)
            temp_out$N2 <- as.factor(N2)
            temp_out$relation <- as.factor("n1_xmap_n2")
            
            pval <- as.numeric((n2v1[[i]])$sl)
            temp_out2 <- data.frame(pval, percent)
            temp_out2$N1 <- as.factor(N1)
            temp_out2$N2 <- as.factor(N2)
            temp_out2$relation <- as.factor("n2_xmap_n1")
            temp_out <- rbind(temp_out, temp_out2)
            
            if (i == 1){
              output <- temp_out
            }
            else{
              output <- rbind(output, temp_out)
            }
            
          }
          
          if (n1 == 1 & n2 == 2){
            output_final <- output
          }
          else{
            output_final <- rbind(output_final, output)
          }
          
          # continously fills up as more variables get cross mapped
          variables.done <- c(variables.done, c=N1)
        }
        
        
      }
      else{
        print("same variable")
      }
    }
  }
  
  # filter out interactions that had no significant (monotonic) trend (by ManKendall)
  # when looking at rho ~ libsize plot (at all quartiles)
  filtered <- filter(output_final, pval <= 0.05)
  
  # now only those significant in all four quartiles are considered
  filtered <- filtered %>% group_by(N1,N2, relation) %>% summarise(percs = length(percent))
  filtered <- filter(filtered, percs > 3) %>% arrange(N1, relation) %>% droplevels()
  
  # now create a list of variables that a focal variable interacts with
  
  xmap <- levels(filtered$relation)
  
  for (i in 1:length(xmap)){
    
    temp <- filter(filtered, relation == xmap[i]) %>% droplevels()
    
    # working with n1_xmap_n2 (so N2 influences N1)
    if (i == 1){
      levs <- levels(temp$N1)
      for (j in 1:length(levs)){
        temp.temp <- filter(temp, N1 == levs[j]) %>% droplevels()
        focal <- levs[j]
        interacting.var <- paste(levels(temp.temp$N2), collapse = "/") 
        direction <- as.factor("focal.influenced.by")
        temp.out <- data.frame(focal, interacting.var, direction)
        if (j == 1){
          ult.out <- temp.out
        }
        else{
          ult.out <- rbind(ult.out, temp.out)
        }
      }
    }
    else {
      levs <- levels(temp$N1)
      for (j in 1:length(levs)){
        temp.temp <- filter(temp, N1 == levs[j]) %>% droplevels()
        focal <- levs[j]
        interacting.var <- paste(levels(temp.temp$N2), collapse = "/") 
        direction <- as.factor("focal.influences")
        temp.out <- data.frame(focal, interacting.var, direction)
        ult.out <- rbind(ult.out, temp.out)
      }
    }
  }

  ult.out$ID <- IDs[ids]
  
  if (ids == 1){
    def.ult.out <- ult.out
  }
  else{
    def.ult.out <- rbind(def.ult.out, ult.out)
  }
  
}
end.time <- Sys.time()

end.time - start.time


setwd("DATA/")

INDIV.CCM <- def.ult.out
save(INDIV.CCM, file = "INDIV.CCM.RData")
variables <- levels(INDIV.CCM$focal)




########################################################
#### PERFORMING CCM FOR STITCHED TS ####################
########################################################
rm(list=ls())

load("STITCHED.TS.TRANSF.RData")

incubators <- levels(STITCHED.TS.TRANSF$incubator)

start.time <- Sys.time()
for (incs in 1:length(incubators)){
  
  
  # filtering out one incubator and make horizontal dfs
  temp <- filter(STITCHED.TS.TRANSF, incubator == incubators[incs]) %>%
    ungroup() %>% 
    select(day, variable, scaled) %>%
    spread(key = variable, value = scaled) %>%
    ungroup() %>%
    select(-day)
  
  # design sequence of library size
  maxlib <- length(temp$temperature)
  libs <- c(round(seq(0.1, 1, 0.1)*maxlib))
  perc <- c("25%","50%","75%","100%")
  
  
  ## CCM ANALYSIS testing N1 cross-mapping N2 (test N2 as causation for N1)
  
  # go through all species combinations
  # generate a df where all variable are crossmapped against each other 
  # loop always checks that no redundant calculations are performed
  
  variables.done <- list(0)
  for (n1 in 1:ncol(temp)){
    for (n2 in 1:ncol(temp)){
      
      N1 <- names(temp)[n1]
      N2 <- names(temp)[n2]
      
      if (n1 != n2){
        
        # ensures that no reduntant calculations are performed
        if (N2 %in% variables.done){
          print("We've done that!")
        }
        else{
          print(paste0("CURRENT INCUBATOR: ", incubators[incs]))
          print(paste0("RUN:", n1))
          print(paste0("Variable N1: ", N1))
          print(paste0("Variable N2: ", N2))
          
          
          # determine embedding dimension
          E.test.n1 = NULL
          n12q=NULL
          for(E.t in 2:8){
            cmxy.t <- ccm(temp, E = E.t, lib_column = N1,
                          target_column = N2, lib_sizes = maxlib,
                          num_samples = 1,
                          tp=-1, random_libs = F, silent = T)
            E.test.n1 = rbind(E.test.n1,cmxy.t)}
          E_n1 <- E.test.n1$E[which.max(E.test.n1$rho)[1]]
          
          # CCM analysis with varying library size (L)
          n1_xmap_n2 <- ccm(temp, E=E_n1,lib_column= N1,
                            target_column= N2,
                            lib_sizes=libs, num_samples=200, replace=T,
                            RNGseed = 2301, silent = T)
          n1_xmap_n2 <- na.omit(n1_xmap_n2)
          n12q= as.matrix(aggregate(n1_xmap_n2[,c('rho')],
                                    by = list(as.factor(n1_xmap_n2$lib_size)),
                                    quantile)[,'x'])
          # ensures that no inf values sneak in (happens only for didinium, diss. oxygen interaction)
          n12q <- n12q[is.finite(rowSums(n12q)),]
          n1v2 <- apply(n12q[,2:5],2,MannKendall)
          
          ##################################s######################
          ## CCM ANALYSIS testing N2 cross-mapping N1 (test N1 as causation for N2)
          E.test.n2 = NULL
          n21q=NULL
          for(E.t in 2:8){
            cmxy.t <- ccm(temp, E = E.t, lib_column = N2,
                          target_column = N1, lib_sizes = maxlib,
                          num_samples = 1,
                          tp=-1, random_libs = F, silent = T)
            E.test.n2 = rbind(E.test.n2,cmxy.t)}
          E_n2 <- E.test.n2$E[which.max(E.test.n2$rho)[1]]
          
          # CCM
          n2_xmap_n1 <- ccm(temp, E=E_n2,lib_column=N2,
                            target_column=N1, 
                            lib_sizes=libs, num_samples=200, replace=T,
                            RNGseed=2301, silent = T)
          n2_xmap_n1 <- na.omit(n2_xmap_n1)
          n21q=as.matrix(aggregate(n2_xmap_n1[,c('rho')],
                                   by = list(as.factor(n2_xmap_n1$lib_size)),
                                   quantile)[,'x'])
          n21q <- n21q[is.finite(rowSums(n21q)),]
          n2v1 <- apply(n21q[,2:5],2,MannKendall)
          
          for (i in 1:4){
            pval <- as.numeric((n1v2[[i]])$sl)
            percent <- as.factor(perc[i])
            temp_out <- data.frame(pval, percent)
            temp_out$N1 <- as.factor(N1)
            temp_out$N2 <- as.factor(N2)
            temp_out$relation <- as.factor("n1_xmap_n2")
            
            pval <- as.numeric((n2v1[[i]])$sl)
            temp_out2 <- data.frame(pval, percent)
            temp_out2$N1 <- as.factor(N1)
            temp_out2$N2 <- as.factor(N2)
            temp_out2$relation <- as.factor("n2_xmap_n1")
            temp_out <- rbind(temp_out, temp_out2)
            
            if (i == 1){
              output <- temp_out
            }
            else{
              output <- rbind(output, temp_out)
            }
            
          }
          
          if (n1 == 1 & n2 == 2){
            output_final <- output
          }
          else{
            output_final <- rbind(output_final, output)
          }
          
          # continously fills up as more variables get cross mapped
          variables.done <- c(variables.done, c=N1)
        }
        
        
      }
      else{
        print("same variable")
      }
    }
  }
  
  # filter out interactions that had no significant (monotonic) trend (by ManKendall)
  # when looking at rho ~ libsize plot (at all quartiles)
  filtered <- filter(output_final, pval <= 0.05)
  
  # now only those significant in all four quartiles are considered
  filtered <- filtered %>% group_by(N1,N2, relation) %>% summarise(percs = length(percent))
  filtered <- filter(filtered, percs > 3) %>% arrange(N1, relation) %>% droplevels()
  
  # now create a list of variables that a focal variable interacts with
  
  xmap <- levels(filtered$relation)
  
  for (i in 1:length(xmap)){
    
    temp <- filter(filtered, relation == xmap[i]) %>% droplevels()
    
    # working with n1_xmap_n2 (so N2 influences N1)
    if (i == 1){
      levs <- levels(temp$N1)
      for (j in 1:length(levs)){
        temp.temp <- filter(temp, N1 == levs[j]) %>% droplevels()
        focal <- levs[j]
        interacting.var <- paste(levels(temp.temp$N2), collapse = "/") 
        direction <- as.factor("focal.influenced.by")
        temp.out <- data.frame(focal, interacting.var, direction)
        if (j == 1){
          ult.out <- temp.out
        }
        else{
          ult.out <- rbind(ult.out, temp.out)
        }
      }
    }
    else {
      levs <- levels(temp$N1)
      for (j in 1:length(levs)){
        temp.temp <- filter(temp, N1 == levs[j]) %>% droplevels()
        focal <- levs[j]
        interacting.var <- paste(levels(temp.temp$N2), collapse = "/") 
        direction <- as.factor("focal.influences")
        temp.out <- data.frame(focal, interacting.var, direction)
        ult.out <- rbind(ult.out, temp.out)
      }
    }
  }
  
  ult.out$incubator <- incubators[incs]
  
  if (incs == 1){
    def.ult.out <- ult.out
  }
  else{
    def.ult.out <- rbind(def.ult.out, ult.out)
  }
  
}
end.time <- Sys.time()

end.time - start.time







for (i in 1:length(incubators)){
  for (j in 1:length(variables)){
  temp <- filter(INDIV.CCM, incubator == incubators[i] & focal == variables[j]) %>% droplevels()
  temp$interacting.var <- as.vector(temp$interacting.var)
  
  focal <- variables[j]
  interacting <- paste0(c(temp$interacting.var), collapse = "/")
  incubator <- incubators[i]
  temp.out <- data.frame(focal, interacting, incubator)
  
  if (i == 1 & j == 1){
    fin.out <- temp.out
  }
  else{
    fin.out <- rbind(fin.out, temp.out)
  }
  
  }
}
end.time <- Sys.time()
end.time - start.time


STITCHED.CCM <- def.ult.out
save(STITCHED.CCM, file = "STITCHED.CCM.RData")
