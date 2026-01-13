# Shannon McWaters, James J. Kearsley, David W. Kikuchi, 
# Timothy J. Polnaszek, Anna Dornhaus
# R script for data analysis and figures for 
# Innovation project: bumble bees, inherent traits, and environment

# Paper reference: XXX (to be updated on acceptance of the paper)
### PAPER VERSION    DRAFT

### Libraries and graphics setup ------------------
library(lme4) #needed for GLMMs
library(lmerTest) #needed for obtaining p-values in lmm
library(emmeans) #post-hoc comparisons
library(ggplot2) #for plots
library(dplyr)

## Colorpalette
simplecomplexcolors <- c("#0072B2", "#D55E00")

####Data Wrangling ---------------------------------------
# Raw data are csv file in same folder
rawdat <- read.csv("BeeInnovationDataRaw.csv")

# Bee's that need to be removed from analysis
## !!!!!!Why
bees_to_remove <- c("A69", "A10")
# Filter the dataframe to remove rows with BeeID values matching to_remove
beeTrimmed <- rawdat %>%
  filter(!BeeID %in% bees_to_remove)

# Main data frame for analysis; same number of rows, limited columns
innovatedata <- beeTrimmed[,c('BeeID','Env','colony')]
#BeeID: unique id number for each bee
#env: whether the bee was subject to a simple or complex environment
innovatedata$Env <- factor(innovatedata$Env, levels = c("s", "c"))
#colony: colony id
no_bees <- length(unique(innovatedata$BeeID))

#Flower1search: time to alight on first novel flower
#Flower1handling: time spent drinking (head in first novel flower)
#resp: responsiveness
innovatedata$resp <- (beeTrimmed$T5travel - beeTrimmed$T4travel)/beeTrimmed$T4travel #differences between travel times for flowers 4 and 5
innovatedata$F1time <- beeTrimmed$Flower1handling + beeTrimmed$Flower1search #total time for 1st novel flower
innovatedata$F1solve <- as.numeric(as.factor(beeTrimmed$Flower1Solve))-1 #binary: success on 1st novel flower or not?
innovatedata$F2time <- beeTrimmed$Flower2handling + beeTrimmed$Flower2search #total time on 2nd novel flower
innovatedata$F2solve <- as.numeric(as.factor(beeTrimmed$Flower2Solve))-1 #binary: success on 2nd novel flower
innovatedata$cap1time <- beeTrimmed$Cap1.handling #handling time on 3rd novel flower
innovatedata$cap1solve <- as.numeric(as.factor(beeTrimmed$Cap1solve))-1 #binary: success on 3rd novel flower
innovatedata$cap2time <- beeTrimmed$Cap2handling #handling time on 4th novel flower
innovatedata$cap2solve <- as.numeric(as.factor(beeTrimmed$Cap2solve))-1 #binary: success on 4th novel flower


# Extract columns for specific flower types
Bumpy <- innovatedata[,c('F1time', 'F1solve','Env','resp')]
names(Bumpy) <- c('time','solve','env','resp')
Folded <- innovatedata[,c('F2time', 'F2solve','Env','resp')]
names(Folded) <- c('time','solve','env','resp')
Cap1 <- innovatedata[,c('cap1time', 'cap1solve','Env','resp')]
names(Cap1) <- c('time','solve','env','resp')
Cap2 <- innovatedata[,c('cap2time', 'cap2solve','Env','resp')]
names(Cap2) <- c('time','solve','env','resp')

# Collage all the data for each flower type in one long data frame, i.e. separate row for each flower type
innovatedatalong <- rbind(Bumpy,Folded,Cap1,Cap2)
innovatedatalong$BeeID <- rep(1:no_bees,4)
innovatedatalong$trial <- rep(c("Bumpy", "Fold", "Cap1", "Cap2"), each = no_bees)
innovatedatalong$trial <- factor(as.factor(innovatedatalong$trial), levels = c("Bumpy", "Fold", "Cap1", "Cap2"))

# removing observations where the bee did not visit the flower
innovatedatalong <- innovatedatalong %>%
  filter(time != 0)

#make data frame for just the time to solve 
innovation <- innovatedatalong %>%
  filter(solve != 0)
#make data frame for just the time to give up
abandoning <- innovatedatalong %>%
  filter(solve != 1)

## ANALYSES -------------------------------------------

####Did the bee solve?####
# env is simple or complex
# trial is flower type
# resp is responsiveness
solve_mod<- glmer(solve ~ env + trial + (1|BeeID), 
                  data = innovatedatalong, 
                  family = binomial)
summary(solve_mod)

# Obtain estimated marginal means (EMMs) for 'trial'
solvemeans <- emmeans(solve_mod, ~ trial)
# Perform pairwise comparisons between trials
pairwise_comparisons <- contrast(solvemeans, method = "pairwise")
# Summary of pairwise comparisons
summary(pairwise_comparisons)

# Result all comparisons n.s.

####Time to solve####
# Descriptive:
hist(innovation$time)
# Note that time to solve has a steeply decreasing, long tail distribution -
# suitable for a log transformation

innovation$logtime <- log(innovation$time)
innovation_mod <- lmer(logtime ~ env + trial + (1|BeeID),
                        data = innovation)
summary(innovation_mod)

# Obtain estimated marginal means (EMMs) for 'trial'
innovmeans <- emmeans(innovation_mod, ~ trial)
# Perform pairwise comparisons between trials
pairwise_comparisons2 <- contrast(innovmeans, method = "pairwise")
# Summary of pairwise comparisons
summary(pairwise_comparisons2)

# Results: 
#  env sign - complex significantly different than simple
#  trial (=flower type) sign
#   bumpy=folded but significantly different from cap1=cap2


####Time to give up####
# Descriptive:
hist(abandoning$time)
# Note that time to give up seems bimodal, with some very short and some very long times

abandoning_mod <- lmer(time ~ env + trial + (1|BeeID), 
                       data = abandoning)
summary(abandoning_mod)
#removing trial as the trials are unevenly represented and there is only one observation for one trial
abandoning_mod2 <- lm(time ~ env + (1|BeeID),
                       data = abandoning)
summary(abandoning_mod2)

# Result: either way env not significant
# trial is, but unevenly represented

