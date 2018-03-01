## PBFE - CONVERGENT CROSS MAPPING OF SPECIES INTERACTIONS
## David Inauen - 01.03.2018
## Mostly based on Chang 2016 - EDM for beginners 
## and supp. material


rm(list=ls())

library(rEDM)
library(Kendall)
library(dplyr)
library(tidyr)


# exploration

load("DATA/INDIV.TS.TRANSF.RData")



# filtering out one ID and make horizontal dfs
temp <- filter(INDIV.TS.TRANSF, ID == "B12") %>%
  ungroup() %>% 
  select(day, variable, differenced) %>%
  spread(key = variable, value = differenced) %>%
  ungroup() %>%
  select(-day)

# design sequence of library size
maxlib <- length(temp$temperature)
libs <- c(round(seq(0.1, 1, 0.1)*maxlib))
perc <- c("25%","50%","75%","100%")


## CCM ANALYSIS testing N1 cross-mapping N2 (test N2 as causation for N1)

# go through all species combinations
for (n1 in 1:ncol(temp)){
  for (n2 in 1:ncol(temp)){
    if (n1 != n2){
      
      N1 <- names(temp)[n1]
      N2 <- names(temp)[n2]
      
      
      # determine embedding dimension
      E.test.n1 = NULL
      n12q=NULL
      for(E.t in 2:8){
        cmxy.t <- ccm(temp, E = E.t, lib_column = N1,
                      target_column = N2, lib_sizes = maxlib,
                      num_samples = 1,
                      tp=-1, random_libs = F)
        E.test.n1 = rbind(E.test.n1,cmxy.t)}
      E_n1 <- E.test.n1$E[which.max(E.test.n1$rho)[1]]
      
      # CCM analysis with varying library size (L)
      n1_xmap_n2 <- ccm(temp, E=E_n1,lib_column= N1,
                        target_column= N2,
                        lib_sizes=libs, num_samples=200, replace=T,
                        RNGseed = 2301)
      n1_xmap_n2 <- na.omit(n1_xmap_n2)
      n12q= as.matrix(aggregate(n1_xmap_n2[,c('rho')],
                                by = list(as.factor(n1_xmap_n2$lib_size)),
                                quantile)[,'x'])
      n1v2 <- apply(n12q[,2:5],2,MannKendall)
      
      ##################################s######################
      ## CCM ANALYSIS testing N2 cross-mapping N1 (test N1 as causation for N2)
      E.test.n2 = NULL
      for(E.t in 2:8){
        cmxy.t <- ccm(temp, E = E.t, lib_column = N2,
                      target_column = N1, lib_sizes = maxlib,
                      num_samples = 1, tp=-1, random_libs = F)
        E.test.n2=rbind(E.test.n2,cmxy.t)}
      E_n2 <- E.test.n2$E[which.max(E.test.n2$rho)[1]]
      
      # CCM
      n2_xmap_n1 <- ccm(temp, E=E_n2,lib_column=N2,
                        target_column=N1, 
                        lib_sizes=libs, num_samples=200, replace=T,
                        RNGseed=2301)
      n2_xmap_n1 <- na.omit(n1_xmap_n2)
      n21q=as.matrix(aggregate(n2_xmap_n1[,c('rho')],
                               by = list(as.factor(n2_xmap_n1$lib_size)),
                               quantile)[,'x'])
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
      
    }
    else{
      print("same variable")
    }

  }
}




filtered <- filter(output_final, pval <= 0.05)

xmaps <- levels(filter$relation)

# now who influences whom:
for(i in 1:2){
  
  temp <- filtered %>% filter(relation == xmaps[i]) %>% droplevels()
  
  if (i == 1){
    
    # first round N2 influences N1 then vice versa!
    levs <- as.factor(levels(temp$N2))
    for (j in 1:length(levs)){
      temp_temp <- temp %>% filter(N2 == levs[j]) %>% droplevels()
      focal <- as.character(levs[j])
      influences <- paste(unique(temp_temp$N1), collapse = "/")
      temp_output <- data.frame(focal, influences)
      if (i == 1 & j == 1){
        output <- temp_output
      }
      else{
        output <- rbind(output, temp_output)
      }
    }

  }
  
  else{
    # second round N1 influences N2 then vice versa!
    levs <- as.factor(levels(temp$N1))
    for (j in 1:length(levs)){
      temp_temp <- temp %>% filter(N1 == levs[j]) %>% droplevels()
      focal <- as.character(levs[j])
      influences <- paste(unique(temp_temp$N2), collapse = "/")
      temp_output <- data.frame(focal, influences)
      output <- rbind(output, temp_output)
    }
  }
  
}




# 
# 
# PLOTTING
# # Plot forecast skill vs library size
# # Plot N1 cross-mapping N2
# plot(n12q[,3]~libs,type="l",col="red",ylim=c(0,1),lwd=2,
#      main="Convergent cross mapping CCM",xlab="Library size",ylab=expression(rho))
# # median predictive skill vs library size (or we can use mean predictive skill)
# lines(n12q[,2]~libs,col="red",lwd=1,lty=2) # 1st quantile
# lines(n12q[,4]~libs,col="red",lwd=1,lty=2) # 3rd quantile
# 
# # Plot N2 cross-mapping N1
# lines(n21q[,3]~libs,col="blue",lwd=1,lty=1) # median
# lines(n21q[,2]~libs,col="blue",lwd=1,lty=2) # 1st quantile
# lines(n21q[,4]~libs,col="blue",lwd=1,lty=2) # 3rd quantile
# legend(10,1,c(paste(N1,N2, sep = " xmap "),
#               paste(N2,N1, sep = " xmap ")),
#        lty=c(1,1),col=c("red","blue"))
# abline(h=cor(temp[,N1],temp[,N2]),lty=3)
# # explanation:
# # red curve gives the causation: second (N2), explains first (N1)
