rawdata <- read.csv(url("https://raw.githubusercontent.com/shannonmcwaters/Innovation/refs/heads/main/BeeInnovationDataRaw.csv"))
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
beeTrimmed$HB10 <- ifelse(!is.na(beeTrimmed$T5_HB), 1, 0)
innovatedata <- beeTrimmed[,c('BeeID','Env','colony')]
innovatedata$F1search <- beeTrimmed$Flower1search
innovatedata$F2search <- beeTrimmed$Flower2search
innovatedata$SRI <- beeTrimmed$SRI
innovatedata$resp <- (beeTrimmed$T5travel - beeTrimmed$AvgPre)/beeTrimmed$AvgPre #differences between travel times for flowers 4 and 5
innovatedata$F1time <- beeTrimmed$Flower1handling   #total time for 1st novel flower
innovatedata$F1solve <- as.numeric(as.factor(beeTrimmed$Flower1Solve))-1 #binary: success on 1st novel flower or not?
innovatedata$F2time <- beeTrimmed$Flower2handling #total time on 2nd novel flower
innovatedata$F2solve <- as.numeric(as.factor(beeTrimmed$Flower2Solve))-1 #binary: success on 2nd novel flower
innovatedata$cap1time <- beeTrimmed$Cap1.handling #handling time on 3rd novel flower
innovatedata$cap1solve <- as.numeric(as.factor(beeTrimmed$Cap1solve))-1 #binary: success on 3rd novel flower
innovatedata$cap2time <- beeTrimmed$Cap2handling #handling time on 4th novel flower
innovatedata$cap2solve <- as.numeric(as.factor(beeTrimmed$Cap2solve))-1 #binary: success on 4th novel flower
innovatedata$HB10 <- beeTrimmed$HB10
innovatedata$H_F1_T1 <- beeTrimmed$H_F1_T1

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
  select(BeeID, trial, time, solve, resp, HB10, SRI, Env, colony, F1search, F2search,H_F1_T1)

# removing observations where the bee did not visit the flower
innovatedatalongA <- innovatedatalong %>%
  filter(time != 0)

# Calculate the proportion of 'no' values for each solve column
prop_not_solved <- colSums(innovatedatalongA[, "solve"] == "0") / nrow(innovatedatalongA)
# Calculate mean time for each trial type
mean_times <- innovatedatalongA %>%
  group_by(trial) %>%
  summarise(mean_time = mean(time, na.rm = TRUE))

#making datafram for those who solved
innovatedatalong_solved <- innovatedatalongA %>%
  filter(solve =="1")
####How does environment play a role####
bee_renamed <- beeTrimmed %>%
  rename(
    Flower1_handling = Flower1handling,
    Flower2_handling = Flower2handling,
    Cap1_handling = Cap1.handling,
    Cap2_handling = Cap2handling,
    Flower1_solve = Flower1Solve,
    Flower2_solve = Flower2Solve,
    Cap1_solve = Cap1solve,
    Cap2_solve = Cap2solve
  )

# Step 2: Pivot to long format
bee_long <- bee_renamed %>%
  select(BeeID, Env,
         Flower1_handling, Flower1_solve,
         Flower2_handling, Flower2_solve,
         Cap1_handling, Cap1_solve,
         Cap2_handling, Cap2_solve) %>%
  pivot_longer(
    cols = -c(BeeID, Env),
    names_to = c("flower", ".value"),
    names_sep = "_"
  ) %>%
  mutate(
    landed = ifelse(handling > 0, 1, 0),
    solve = ifelse(solve == "yes", 1,
                   ifelse(solve == "no", 0, NA))  # clean NAs if present
  )

glm_landing_fixed <- glm(landed ~ Env + flower, family = "binomial", data = bee_long)
summary(glm_landing_fixed)

#env and solving
bee_landed <- bee_long %>% filter(landed == 1)
env_mod = lm(handling~Env+flower,data=bee_landed)
summary(env_mod)


# Logistic regression
glm_solve_all <- glm(solve ~ Env+flower, data = bee_landed, family = "binomial")
summary(glm_solve_all)


####personality and solving####
trait_summary <- innovatedatalong %>%
  mutate(landed = ifelse(time > 0, 1, 0)) %>%
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
    Env = first(Env)
  )

#predcit landing
model_SRI_landed <- lm(prop_landed ~ SRI, data = trait_summary)
model_HB10_landed <- lm(prop_landed ~ HB10, data = trait_summary)
model_resp_landed <- lm(prop_landed ~ resp, data = trait_summary)

summary(model_SRI_landed)
summary(model_HB10_landed)
summary(model_resp_landed)

trait_summaryA <- innovatedatalongA %>%
  mutate(landed = ifelse(time > 0, 1, 0)) %>%
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
    Env = first(Env)
  )
#predict solving
model_SRI_solved <- lm(prop_solved ~ SRI, data = trait_summaryA)
model_HB10_solved <- lm(prop_solved ~ HB10, data = trait_summaryA)
model_resp_solved <- lm(prop_solved ~ resp, data = trait_summaryA)

summary(model_SRI_solved)
summary(model_HB10_solved)
summary(model_resp_solved)

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


####individual traits predicting Innovativeness#### 
#SRI
sri_lm=lm(time~SRI * trial, data=innovatedatalong_solved)
summary(sri_lm)
#Resp
resp_lm = lm(time~resp*trial, data=innovatedatalong_solved)
summary(resp_lm)
#exploration
exp_lm <- lm(time~HB10 * trial,data=innovatedatalong_solved)
summary(exp_lm)

#handling time of first flower
hand_lm <- lm(time~ H_F1_T1 +trial, data=innovatedatalongA)
summary(hand_lm)

####Search time###
bumpy_folded_search <- innovatedatalongA %>%
  filter(trial %in% c("Bumpy", "Folded")) %>%
  mutate(
    search_time = case_when(
      trial == "Bumpy"  ~ F1search,
      trial == "Folded" ~ F2search
    )
  ) %>%
  select(BeeID, trial, search_time, everything())

searchmod <- lm(time ~ search_time + trial, data = bumpy_folded_search)
summary(searchmod)
bee_avg_search <- bumpy_folded_search %>%
  group_by(BeeID, Env) %>%
  summarise(mean_search = mean(search_time, na.rm = TRUE), .groups = "drop")
wilcox.test(mean_search ~ Env, data = bee_avg_search)
####trial anova/ tukey####
trial <- aov(time ~ trial, data = innovatedatalong)
TukeyHSD(trial)





