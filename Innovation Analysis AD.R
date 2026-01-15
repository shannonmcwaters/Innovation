# Shannon McWaters, James J. Kearsley, David W. Kikuchi, 
# Timothy J. Polnaszek, Anna Dornhaus
# R script for data analysis and figures for 
# Innovation project: bumble bees, inherent traits, and environment

# Paper reference: XXX (to be updated on acceptance of the paper)
### PAPER VERSION    DRAFT
### PART I: DATA IMPORT, FORMATTING, AND ANALYSIS

### Libraries ------------------
library(lme4) #needed for GLMMs
library(lmerTest) #needed for obtaining p-values in lmm
library(emmeans) #post-hoc comparisons
library(ggplot2) #for plots
library(dplyr)
library(coxme)
library(stringdist)
library(tidyr)

#### IMPORT DATA -----------------------------------
# Raw data are csv file in same folder
rawdata <- read.csv("BeeInnovationDataRaw.csv")

#### FUNCTION DEFINITIONS --------------------------

# SRI (sequence repetition index for looking at routine-ness,
# specifically how much bees revisit flowers in the exact same order)
calculate_repetition_index <- function(orders) {
  # Remove missing or blank trials
  valid_orders <- orders[!is.na(orders) & orders != ""]
  # If fewer than 2 valid sequences, return NA
  if (length(valid_orders) < 2) return(NA)
  # Calculate pairwise Levenshtein distances
  distance_matrix <- stringdistmatrix(valid_orders, valid_orders, method = "lv")
  # Normalize distances (divide by sequence length)
  sequence_length <- nchar(valid_orders[1])  # Assume all sequences are the same length
  normalized_distances <- distance_matrix / sequence_length
  # Convert distances to similarity scores (1 - distance)
  similarity_matrix <- 1 - normalized_distances
  # Calculate the average similarity (excluding diagonal values)
  repetition_index <- mean(similarity_matrix[lower.tri(similarity_matrix)])
  return(repetition_index)
}

#### Data Wrangling ---------------------------------

beedata <- rawdata
# Calculate average time spent traveling over trips 2-4 in the pre-training 
beedata$AvgPre <- rowMeans(beedata[,c("T2_T","T3_T","T4_T")], na.rm=TRUE)

# Calculate SRI for each bee and add it to the temp dataset
beedata <- beedata %>%
  group_by(BeeID) %>%
  mutate(SRI = calculate_repetition_index(c(T1_order, T2_order, T3_order, T4_order)))

# Define exploration as whether or not bees landed on blue flowers in trial 5 
beedata$HB10 <- ifelse(!is.na(beedata$T5_HB), 1, 0)

#BeeID: unique id number for each bee
#colony: colony id
no_bees <- length(unique(beedata$BeeID))

#env: whether the bee was subject to a simple or complex environment
beedata$env <- factor(beedata$Env, levels = c("s", "c"))

# Calculating additional summary statistics
#  Flower1search: time to alight on first novel flower
#  Flower1handling: time spent drinking (head in first novel flower)
#  resp: responsiveness
beedata$F1search <- beedata$Flower1search
beedata$F2search <- beedata$Flower2search
beedata$resp <- (beedata$T5travel - beedata$T4travel)/beedata$T4travel #differences between travel times for flowers 4 and 5
beedata$F1time <- beedata$Flower1handling + beedata$Flower1search #total time for 1st novel flower
beedata$F1solve <- as.numeric(as.factor(beedata$Flower1Solve))-1 #binary: success on 1st novel flower or not?
beedata$F2time <- beedata$Flower2handling + beedata$Flower2search #total time on 2nd novel flower
beedata$F2solve <- as.numeric(as.factor(beedata$Flower2Solve))-1 #binary: success on 2nd novel flower
beedata$cap1time <- beedata$Cap1.handling #handling time on 3rd novel flower
beedata$cap1solve <- as.numeric(as.factor(beedata$Cap1solve))-1 #binary: success on 3rd novel flower
beedata$cap2time <- beedata$Cap2handling #handling time on 4th novel flower
beedata$cap2solve <- as.numeric(as.factor(beedata$Cap2solve))-1 #binary: success on 4th novel flower

# Extract columns for specific flower types
Bumpy <- beedata[,c('F1time', 'F1solve','env','resp')]
names(Bumpy) <- c('time','solve','env','resp')
Folded <- beedata[,c('F2time', 'F2solve','env','resp')]
names(Folded) <- c('time','solve','env','resp')
Cap1 <- beedata[,c('cap1time', 'cap1solve','env','resp')]
names(Cap1) <- c('time','solve','env','resp')
Cap2 <- beedata[,c('cap2time', 'cap2solve','env','resp')]
names(Cap2) <- c('time','solve','env','resp')

# Now we construct the main data frame for analysis:
# Collage all the data for each flower type in one long data frame, i.e. separate row for each flower type and bee trip
innovationlong <- beedata %>%
  pivot_longer(cols = c(F1time, F2time, cap1time, cap2time), 
               names_to = "trial", 
               values_to = "time") %>%
  mutate(trial = recode(trial, 
                        "F1time" = "Bumpy", 
                        "F2time" = "Folded", 
                        "cap1time" = "Cap1", 
                        "cap2time" = "Cap2")) %>%
  # Add solve column using corresponding solve columns
  mutate(solve = case_when(
    trial == "Bumpy" ~ F1solve,
    trial == "Folded" ~ F2solve,
    trial == "Cap1" ~ cap1solve,
    trial == "Cap2" ~ cap2solve
  )) %>%
  # Keep only necessary columns
  select(BeeID, trial, time, solve, resp, HB10, SRI, env, colony, F1search, F2search,H_F1_T1)

innovationlong$trial <- factor(innovationlong$trial, levels = c("Bumpy", "Folded", "Cap1", "Cap2"))

# Descriptive:
hist(innovationlong$time)
# Note that time to solve has a steeply decreasing, long tail distribution -
# suitable for a log transformation
innovationlong$logtime <- log(innovationlong$time)

# Removing observations where the bee did not visit the flower
innovationlong$landed <- ifelse(innovationlong$time != 0, TRUE, FALSE)
innovationlanded <- innovationlong %>%
  filter(landed)

# Make data frame for just the time to solve 
innovationsuccess <- innovationlanded %>%
  filter(solve != 0)

abandonedinnovation <- innovationlanded %>%
  filter(solve != 1)

#### SUMMARIES ------------------------------

# Calculate the proportion of 'no' values for each solve column
prop_not_solved <- colSums(innovationlanded[, "solve"] == "0") / nrow(innovationlanded)

# Calculate mean solving time for each trial type
mean_times <- innovationsuccess %>%
  group_by(trial) %>%
  summarise(mean_time = mean(time, na.rm = TRUE))

# Search time in first two innovation trials
bumpy_folded_search <- innovationlanded %>%
  filter(trial %in% c("Bumpy", "Folded")) %>%
  mutate(
    search_time = case_when(
      trial == "Bumpy"  ~ F1search,
      trial == "Folded" ~ F2search
    )
  ) %>%
  select(BeeID, trial, search_time, everything())

## ANALYSES -------------------------------------------

#### EFFECT OF FLOWER TYPE ----------------------------
trialmod <- aov(logtime ~ trial, data = innovationsuccess)
summary(trialmod)
TukeyHSD(trialmod)


#### EFFECT OF ENVIRONMENT ----------------------------

# On landing:
glm_landing_fixed <- glm(landed ~ env + trial, family = "binomial", data = innovationlong)
summary(glm_landing_fixed)
# Since there are only 3 cases of no landing, this is probably not a useful analysis.

# On solving:
glm_solve_all <- glm(solve ~ env + trial, data = innovationlanded, family = "binomial")
summary(glm_solve_all)

# On time to solve:
env_mod <- lm(logtime ~ env + trial, data = innovationlanded)
summary(env_mod)

# On time to give up:
# Descriptive:
hist(abandonedinnovation$time)
# Note that time to give up seems bimodal, with some very short and some very long times
abandoning_mod <- lm(time ~ env + trial, 
                       data = abandonedinnovation)
summary(abandoning_mod)

# On search time:
bee_avg_search <- bumpy_folded_search %>%
  group_by(BeeID, env) %>%
  summarise(mean_search = mean(search_time, na.rm = TRUE), .groups = "drop")

wilcox.test(mean_search ~ env, data = bee_avg_search)


#### EFFECT OF BEE TRAITS ----------------------------

# For landing (yes/no) and solving (yes/no), we use one data point per bee, to 
# calculate proportion of time landing/solving. 

# Make summary data table where each bee is only one row again
trait_summary <- innovationlong %>%
  group_by(BeeID) %>%
  summarise(
    prop_landed = mean(landed, na.rm = TRUE),
    prop_solved = if (sum(landed) > 0) {
      mean(solve[landed == 1], na.rm = TRUE)
    } else {
      NA_real_
    },
    SRI = first(SRI),
    HB10 = first(HB10),
    resp = first(resp),
    env = first(env)
  )

# On landing:
model_SRI_landed <- lm(prop_landed ~ SRI, data = trait_summary)
model_HB10_landed <- lm(prop_landed ~ HB10, data = trait_summary)
model_resp_landed <- lm(prop_landed ~ resp, data = trait_summary)
summary(model_SRI_landed)
summary(model_HB10_landed)
summary(model_resp_landed)
# Again, not a lot of cases of not landing, so not really expecting an effect

# On solving:
model_SRI_solved <- lm(prop_solved ~ SRI, data = trait_summary)
model_HB10_solved <- lm(prop_solved ~ HB10, data = trait_summary)
model_resp_solved <- lm(prop_solved ~ resp, data = trait_summary)
summary(model_SRI_solved)
summary(model_HB10_solved)
summary(model_resp_solved)

# On time to solve: 
# For the continuous outcome of solving time, we use all innovation trials 
# (flower types) separately, thus multiple measurements per bee; 
# but only trials in which the bee actually solved (i.e. accessed reward).

# SRI (routine formation)
sri_lm <- lm(logtime ~ SRI + trial, data = innovationsuccess)
summary(sri_lm)

# Responsiveness
resp_lm <- lm(logtime ~ resp + trial, data = innovationsuccess)
summary(resp_lm)

# Exploration
exp_lm <- lm(logtime ~ HB10 + trial,data = innovationsuccess)
summary(exp_lm)

# Handling time of first flower
hand_lm <- lm(logtime ~ log(H_F1_T1) + trial, data = innovationsuccess)
summary(hand_lm)

# Search time in first two innovation trials on solving time in those trials
searchmod <- lm(logtime ~ log(search_time) + trial, data = bumpy_folded_search)
summary(searchmod)

# Individual ID on solving time - i.e. some indivdiual traits not captured above
BeeID_mod <- lmer(logtime ~ trial + (1 | BeeID), data = innovationsuccess)
BeeID_mod2 <- lm(logtime ~ trial, data = innovationsuccess)
anova(BeeID_mod, BeeID_mod2)

# Get variance components
var_components <- VarCorr(BeeID_mod)
# Extract the variance components
var_beeID <- var_components$BeeID[1]
var_residual <- attr(var_components, "sc")^2
# Calculate repeatability
repeatability <- var_beeID / (var_beeID + var_residual)
repeatability


#### ALT Analysis ----------------------------------------

## Putting it all in together does not change result
# On solving:
solve_mod<- glmer(solve ~ env + trial + (1|BeeID), 
                  data = innovationlanded, 
                  family = binomial)
summary(solve_mod)
# Obtain estimated marginal means (EMMs) for 'trial'
solvemeans <- emmeans(solve_mod, ~ trial)
# Perform pairwise comparisons between trials
pairwise_comparisons <- contrast(solvemeans, method = "pairwise")
# Summary of pairwise comparisons
summary(pairwise_comparisons)
# Result all comparisons n.s.

# On time to solve:
innovation_mod <- lmer(logtime ~ env + trial + (1|BeeID),
                        data = innovationsuccess)
summary(innovation_mod)
# Obtain estimated marginal means (EMMs) for 'trial'
innovmeans <- emmeans(innovation_mod, ~ trial)
# Perform pairwise comparisons between trials
pairwise_comparisons2 <- contrast(innovmeans, method = "pairwise")
# Summary of pairwise comparisons
summary(pairwise_comparisons2)

