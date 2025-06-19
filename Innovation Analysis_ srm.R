rawdata <- read.csv(url("https://raw.githubusercontent.com/shannonmcwaters/Innovation/refs/heads/main/BeeInnovationData_full.csv"))
####packages####
library(lme4) #needed for GLMMs
library(lmerTest) #needed for obtaining p-values in lmm
library(emmeans) #post-hoc comparisons
library(ggplot2)#for plots
library(dplyr)
library(coxme)
library(stringdist)
library(tidyr)
####Data Wrangling####
#Make a new dataframe for wrangling
beeTrimmed <- rawdata
#calculate average time spent traveling over trips 2-4 in the pre-training 
beeTrimmed$AvgPre <- rowMeans(beeTrimmed[,c("T2_T","T3_T","T4_T")], na.rm=TRUE)
#SRI (sequence repetition index for looking at reoutine-ness)
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
# Calculate SRI for each bee and add it to the temp dataset
beeTrimmed <- beeTrimmed %>%
  group_by(BeeID) %>%
  mutate(SRI = calculate_repetition_index(c(T1_order, T2_order, T3_order, T4_order)))

innovatedata <- beeTrimmed[,c('BeeID','Env','colony')]
innovatedata$F1search <- beeTrimmed$Flower1search
innovatedata$F2search <- beeTrimmed$Flower2search
innovatedata$SRI <- beeTrimmed$SRI
innovatedata$resp <- (beeTrimmed$T5_T - beeTrimmed$AvgPre)/beeTrimmed$AvgPre #differences between travel times for flowers 4 and 5
innovatedata$F1time <- beeTrimmed$Flower1handling   #total time for 1st novel flower
innovatedata$F1solve <- as.numeric(as.factor(beeTrimmed$Flower1Solve))-1 #binary: success on 1st novel flower or not?
innovatedata$F2time <- beeTrimmed$Flower2handling #total time on 2nd novel flower
innovatedata$F2solve <- as.numeric(as.factor(beeTrimmed$Flower2Solve))-1 #binary: success on 2nd novel flower
innovatedata$cap1time <- beeTrimmed$Cap1.handling #handling time on 3rd novel flower
innovatedata$cap1solve <- as.numeric(as.factor(beeTrimmed$Cap1solve))-1 #binary: success on 3rd novel flower
innovatedata$cap2time <- beeTrimmed$Cap2handling #handling time on 4th novel flower
innovatedata$cap2solve <- as.numeric(as.factor(beeTrimmed$Cap2solve))-1 #binary: success on 4th novel flower
innovatedata$HB10 <- ifelse(!is.na(beeTrimmed$T5_HB), 1, 0)

innovatedatalong <- innovatedata %>%
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
  select(BeeID, trial, time, solve, resp, HB10, SRI, Env, colony, F1search, F2search)

# removing observations where the bee did not visit the flower
innovatedatalong <- innovatedatalong %>%
  filter(time != 0)

# Calculate the proportion of 'no' values for each solve column
prop_not_solved <- colSums(innovatedatalong[, "solve"] == "0") / nrow(innovatedatalong)
# Calculate mean time for each trial type
mean_times <- innovatedatalong %>%
  group_by(trial) %>%
  summarise(mean_time = mean(time, na.rm = TRUE))

#making datafram for those who solved
innovatedatalong_solved <- innovatedatalong %>%
  filter(solve =="1")
####How does environment play a role in innovation####
envmod = lm(time ~ Env + trial, data=innovatedatalong)
summary(envmod)
####Innovativeness as a personality trait####
# Fit model with BeeID as random effect
BeeID_mod <- lmer(time ~ trial + (1 | BeeID), data = innovatedatalong_solved)
BeeID_mod2 <- lm(time ~ trial, data = innovatedatalong_solved)
anova(BeeID_mod,BeeID_mod2)
# Get variance components
var_components <- VarCorr(BeeID_mod)
# Extract the variance components
var_beeID <- var_components$BeeID[1]
var_residual <- attr(var_components, "sc")^2
# Calculate repeatability
repeatability <- var_beeID / (var_beeID + var_residual)
repeatability

library(rptR)

rpt_result <- rpt(time ~ (1 | BeeID), grname = "BeeID",
                  data = innovatedatalong, datatype = "Gaussian",
                  nboot = 1000, npermut = 1000)

summary(rpt_result)
####individual traits predicting Innovativeness#### 
#SRI
sri_lm=lmer(time~SRI+trial+(1|BeeID), data=innovatedatalong_solved)
summary(sri_lm)
#Resp
resp_lm = lmer(time~resp+trial+(1|BeeID), data=innovatedatalong_solved)
summary(resp_lm)
#exploration
exp_lm <- lmer(time~HB10+trial+(1|BeeID),data=innovatedatalong_solved)
summary(exp_lm)

#summary(lm(avgInn~T5_HB + SRI,data=beeTrimmed))

####trial anova/ tukey####
trial <- aov(time ~ trial, data = innovatedatalong)
TukeyHSD(trial)

####Search time####
bee_search <- beeTrimmed %>%
  select(BeeID, Env, 
         Flower1search, Flower2search, 
         Flower1handling, Flower2handling) %>%
  pivot_longer(
    cols = c(Flower1search, Flower2search, Flower1handling, Flower2handling),
    names_to = c("flower", ".value"),
    names_pattern = "(Flower[12])(search|handling)"
  ) %>%
  mutate(
    flower = case_when(
      flower == "Flower1" ~ "Flower 1",
      flower == "Flower2" ~ "Flower 2"
    ),
    total_time = search + handling
  )
search_mod <- lm(search ~ Env + flower, data = bee_search)
summary(search_mod)
search_mod <- lm(total_time ~ Env + flower, data = bee_search)
summary(search_mod)






# Step 1: Sort within bee by trial (assuming trial order is chronological)
handling_pairs <- innovatedatalong %>%
  arrange(BeeID, trial) %>%
  group_by(BeeID) %>%
  slice_head(n = 2) %>%  # get first two trials per bee
  mutate(trial_number = row_number()) %>%
  ungroup()

# Step 2: Pivot to wide format (1 row per bee)
handling_wide <- handling_pairs %>%
  select(BeeID, trial_number, time) %>%
  pivot_wider(names_from = trial_number, values_from = time, names_prefix = "time")

# Spearman correlation (rank-based, more robust)
cor.test(handling_wide$time1, handling_wide$time2, method = "spearman")
cor.test(handling_wide$time3, handling_wide$time4, method = "spearman")



