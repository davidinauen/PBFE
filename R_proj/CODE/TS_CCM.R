## PBFE - CONVERGENT CROSS MAPPING OF SPECIES INTERACTIONS
## David Inauen - 01.03.2018
## Mostly based on Chang 2016 - EDM for beginners 
## and supp. material


rm(list=ls())

library(rEDM)
library(Kendall)
library(dplyr)
library(tidyr)
 



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
