## DATA PREPARATION ##
## David Inauen - 1.03.2018
## Script to prepare Data from the PBFE experiment
## for forecast analysis by EDM


#### PREPARING DATA INTO INDIVIDUAL AND STITCHED TIME SERIES -----
rm(list=ls())

load("RAW.RData")

library(dplyr)
library(haven)
library(pracma)
library(dplyr)
library(tidyr)
library(stats)
library(ggplot2)


# throw out Spirostomum sp. (breaks EDM loops later on, because of very low detectability)
temp.df <- droplevels(filter(pbfe_raw, variable != "spirostomum_sp"))

# filter particles to lump together
temp.df.lump <- filter(pbfe_raw, variable == "green-white_particle" | variable == "small-white_particle")

# lump particles,
# green-white and small-white particles exhibit very similar abundance dynamics, highly correlated
# thus taken together as "small (cryptic) flagellate particles"
temp.df.lump <- temp.df.lump %>% ungroup() %>% 
  group_by(ID, date, treat, incubator, type, day, unit) %>%
  summarise(response = sum(response))

# bind data frames
temp.df.lump$variable <- as.factor("small_particle")
temp.df <- droplevels(filter(temp.df, variable != "green-white_particle" & variable != "small-white_particle"))
temp.df.lump <- temp.df.lump[c("variable", "ID", "date", "treat", "incubator", "response", "type", "day", "unit")]
temp.df.lump <- temp.df.lump %>% filter(variable == "small_particle")
temp.df.lump <- as.data.frame(temp.df.lump)
temp.df.lumped <- rbind(temp.df, temp.df.lump)

# remove incomplete carbon and nitrogen data
temp.df.lumped <- temp.df.lumped %>% 
  filter(variable != "carbon_cent" &
           variable != "carbon_non-cent" &
           variable != "nitrogen_cent" &
           variable != "nitrogen_non-cent") %>%
  droplevels()
  
# some incubators were incorrectly labelled
ID <- temp.df.lumped$ID 
inc_182 <- c("B01", "B02", "B03")
inc_181 <- c("B04", "B05", "B06")
inc_180 <- c("B07", "B08", "B09")
inc_179 <- c("B10", "B11", "B12")
inc_178 <- c("B13", "B14", "B15")
inc_177 <- c("B16", "B17", "B18")
temp.df.lumped$incubator <- ifelse(ID %in% inc_182, "182", ifelse(ID %in% inc_181, "181", ifelse(ID %in% inc_180, "180", ifelse(ID %in% inc_179, "179", ifelse(ID %in% inc_178, "178", "177")))))
temp.df.lumped$incubator <- as.factor(temp.df.lumped$incubator)

# save data frame for individual time series
INDIV.TS <- temp.df.lumped

IDs <- levels(INDIV.TS$ID)
incubators <- levels(INDIV.TS$incubator)
variables <- levels(INDIV.TS$variable)

# create stitched time series using replicates within an incubator 
# time series are stitched after oneanother
# Hsieh 2008 - Extending Nonlinear Analysis to Short Ecological Time Series
# "composing it together should not introduce non-linearity!"

for (i in 1:length(incubators)){
  for (j in 1:length(variables)){
    
    temp <- filter(INDIV.TS, incubator == incubators[i] & variable == variables[j])
    
    temp <- temp %>% arrange(ID, day)
    daylist.ori <- unique(temp$day)
    daylist <- append(daylist.ori, daylist.ori+155)
    daylist <- append(daylist, daylist.ori+310)
    temp$day <- daylist
    
    if (i == 1 & j == 1){
      output <- temp
    }
    else{
      output <- rbind(output, temp)
    }
    
  }
}

STITCHED.TS <- output



#### TRANSFORMING INDIVIDUAL.TS -----
# for now, performing the following transformations (mostly based on: "EDM for beginners" - Chang 2017)
# esp. the supplemental .docx file
# - interpolation (cubic hermite)
# - 4th root (takes care of spiky data)
# - first order differencing ()
# - normalising/standartising

# open for discussion:
# - apply detrend (smoothers are discouraged though, see Chang 2017) instead first order differencing?


n.data.points <- length(unique(INDIV.TS$day))
xout <- seq(min(INDIV.TS$day), max(INDIV.TS$day), length = n.data.points)

# interpolation function (cubic)
Interpolate <- function(x, y) {
  xout.inner <- xout[xout>=min(x) & xout<=max(x)]
  res <- pracma::interp1(x=as.numeric(x),
                         y=as.numeric(y),
                         xi=xout.inner,
                         method="cubic")
  res.df <- data.frame(day=xout.inner,
                       response=res)
}

# sorting
INDIV.TS <- arrange(INDIV.TS, ID, variable, day)
# interpolating
res <- group_by(INDIV.TS, variable, ID) %>%
  do(dummy=Interpolate(.$day, .$response))
INDIV.TS.TRANSF <- unnest(res)

# adding the fourth root transform
INDIV.TS.TRANSF$fr.value <- INDIV.TS.TRANSF$response^0.25

# first order differencing
IDs <- as.factor(levels(INDIV.TS.TRANSF$ID))
variables <- as.factor(levels(INDIV.TS.TRANSF$variable))

for (i in 1:length(IDs)){
  for (j in 1:length(variables)){
    
    temp.df <- filter(INDIV.TS.TRANSF, ID == IDs[i] & variable == variables[j])
    
    differenced <- diff(temp.df$fr.value)
    # making time series one time point smaller
    temp.df <- temp.df[-1,]
    temp.df$differenced <- differenced
    
    if (i == 1 & j == 1){
      output <- temp.df
    }
    else{
      output <- rbind(output, temp.df)
    }
    
  }
}

# rescaling 
INDIV.TS.TRANSF <- output
INDIV.TS.TRANSF <- group_by(INDIV.TS.TRANSF, ID, variable) %>%
  mutate(scaled=as.numeric(scale(differenced)))

# the effect of transforming:
effect <- INDIV.TS.TRANSF %>%
  filter(ID == "B12" & variable == "bacteria") %>%
  gather(key = mode, value = response, response, fr.value,
         differenced, scaled)
effect$mode <- as.factor(effect$mode)
levels(effect$mode)
effect <- as.data.frame(effect)
levels(effect$mode) <- c("C - DIFFERENCED", "B - Fourth root", "A - INTERPOLATED", "D - SCALED")
  
# ugh, ggplot does no longer list things alphabetically?
ggplot(effect, aes(day, response)) +
  geom_line() +
  facet_wrap(~mode, scales = "free_y") +
  theme_minimal() +
  ggtitle("Effect of Transformation (Individual Time series)")

#### TRANSFORMING STITCHED.TS -----


n.data.points <- length(unique(STITCHED.TS$day))
xout <- seq(min(STITCHED.TS$day), max(STITCHED.TS$day), length = n.data.points)

# interpolating
res <- group_by(STITCHED.TS, variable, incubator) %>%
  do(dummy=Interpolate(.$day, .$response))
STITCHED.TS.TRANSF <- unnest(res)

# adding the fourth root transform
STITCHED.TS.TRANSF$fr.value <- STITCHED.TS.TRANSF$response^0.25

# first order differencing
incubators <- as.factor(levels(STITCHED.TS.TRANSF$incubator))
variables <- as.factor(levels(STITCHED.TS.TRANSF$variable))

for (i in 1:length(incubators)){
  for (j in 1:length(variables)){
    
    temp.df <- filter(STITCHED.TS.TRANSF, incubator == incubators[i] & variable == variables[j])
    
    differenced <- diff(temp.df$fr.value)
    # making time series one time point smaller
    temp.df <- temp.df[-1,]
    temp.df$differenced <- differenced
    
    if (i == 1 & j == 1){
      output <- temp.df
    }
    else{
      output <- rbind(output, temp.df)
    }
    
  }
}


# rescaling 
STITCHED.TS.TRANSF <- output
STITCHED.TS.TRANSF <- group_by(STITCHED.TS.TRANSF, incubator, variable) %>%
  mutate(scaled=as.numeric(scale(differenced)))

# the effect of transforming:
effect <- STITCHED.TS.TRANSF %>%
  filter(incubator == "177" & variable == "bacteria") %>%
  gather(key = mode, value = response, response, fr.value,
         differenced, scaled)
effect$mode <- as.factor(effect$mode)
levels(effect$mode)
effect <- as.data.frame(effect)
levels(effect$mode) <- c("C - DIFFERENCED", "B - Fourth root", "A - INTERPOLATED", "D - SCALED")

# ugh, ggplot does no longer list things alphabetically?
ggplot(effect, aes(day, response)) +
  geom_line() +
  facet_wrap(~mode, scales = "free_y") +
  theme_minimal() + 
  ggtitle("Effect of Transformation (Stitched Times series)")

save(INDIV.TS.TRANSF, file = "INDIV.TS.TRANSF.RData")
save(STITCHED.TS.TRANSF, file = "STITCHED.TS.TRANSF.RData")
