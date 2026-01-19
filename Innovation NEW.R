# Shannon McWaters, James J. Kearsley, David W. Kikuchi, 
# Timothy J. Polnaszek, Anna Dornhaus
# R script for data analysis and figures for 
# Innovation project: bumble bees, inherent traits, and environment

# Paper reference: XXX (to be updated on acceptance of the paper)
### PAPER VERSION    DRAFT
### PART I: DATA IMPORT, FORMATTING, AND ANALYSES

### Libraries ------------------
library(lme4) #needed for GLMMs
library(lmerTest) #needed for obtaining p-values in lmm
library(emmeans) #post-hoc comparisons
library(coxme)
library(stringdist)
library(tidyverse)
library(scales)  # for number_format

# Graphics setup -----------------------------

## Colorpalette
simplecomplexcolors <- c("#c7d6d5", "#d4c1e3")
simcomcolors_dark <- c("#5a9a8f", "#7b6ea8")


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
#BeeID: unique id number for each bee
#colony: colony id
no_bees <- length(unique(beedata$BeeID))

# Calculate average time spent traveling over trips 2-4 in the pre-training 
beedata$AvgPre <- rowMeans(beedata[,c("T2_T","T3_T","T4_T")], na.rm=TRUE)

# Calculate SRI for each bee and add it to the temp dataset
beedata <- beedata %>%
  group_by(BeeID) %>%
  mutate(SRI = calculate_repetition_index(c(T1_order, T2_order, T3_order, T4_order)))

# Define exploration as whether or not bees landed on blue flowers in trial 5 
beedata$HB10 <- ifelse(!is.na(beedata$T5_HB), 1, 0)

#env: whether the bee was subject to a simple or complex environment
beedata$env <- factor(beedata$Env, levels = c("s", "c"))

# Calculating additional summary statistics
#  Flower1search: time to alight on first novel flower
#  Flower1handling: time spent drinking (head in first novel flower)
#  resp: responsiveness
beedata$resp <- (beedata$T5travel - beedata$AvgPre)/beedata$AvgPre #differences between travel times for trip 5 and avg of trips before
beedata$F1search <- beedata$Flower1search
beedata$F2search <- beedata$Flower2search
beedata$F1solve <- as.numeric(as.factor(beedata$Flower1Solve))-1 #binary: success on 1st novel flower or not?
beedata$F2solve <- as.numeric(as.factor(beedata$Flower2Solve))-1 #binary: success on 2nd novel flower
beedata$cap1solve <- as.numeric(as.factor(beedata$Cap1solve))-1 #binary: success on 3rd novel flower
beedata$cap2solve <- as.numeric(as.factor(beedata$Cap2solve))-1 #binary: success on 4th novel flower
beedata$F1time <- ifelse(beedata$F1solve == 1, beedata$Flower1handling, NA) # solving time for 1st novel flower
beedata$F2time <- ifelse(beedata$F2solve == 1, beedata$Flower2handling, NA) # solving time on 2nd novel flower
beedata$cap1time <- ifelse(beedata$cap1solve == 1, beedata$Cap1.handling, NA) # solving time on 3rd novel flower
beedata$cap2time <- ifelse(beedata$cap2solve == 1, beedata$Cap2handling, NA) # solving time on 4th novel flower
beedata$landed1 <- ifelse(beedata$Flower1handling>0, 1, 0)
beedata$landed2 <- ifelse(beedata$Flower2handling>0, 1, 0)
beedata$landed3 <- ifelse(beedata$Cap1.handling>0, 1, 0)
beedata$landed4 <- ifelse(beedata$Cap2handling>0, 1, 0)
beedata$F1givinguptime <- ifelse(beedata$landed1 & beedata$F1solve == 0, beedata$Flower1handling, NA) # giving up time for 1st novel flower
beedata$F2givinguptime <- ifelse(beedata$landed2 & beedata$F2solve == 0, beedata$Flower2handling, NA) # giving up time on 2nd novel flower
beedata$cap1givinguptime <- ifelse(beedata$landed3 & beedata$cap1solve == 0, beedata$Cap1.handling, NA) # giving up time on 3rd novel flower
beedata$cap2givinguptime <- ifelse(beedata$landed4 & beedata$cap2solve == 0, beedata$Cap2handling, NA) # giving up time on 4th novel flower

# Calculate the proportion of successful landings and successful solving per bee
beedata$prop_landed <- rowSums(cbind(beedata$landed1, beedata$landed2, beedata$landed3, beedata$landed4), na.rm = TRUE) / 4
beedata$prop_solved <- rowSums(cbind(beedata$F1solve, beedata$F2solve, beedata$cap1solve, beedata$cap2solve), na.rm = TRUE) / 4
beedata$avg_time <- rowMeans(cbind(beedata$F1time, beedata$F2time, beedata$cap1time, beedata$cap2time), na.rm = TRUE)

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
  mutate(landed = case_when(
    trial == "Bumpy" ~ landed1,
    trial == "Folded" ~ landed2,
    trial == "Cap1" ~ landed3,
    trial == "Cap2" ~ landed4
  ))

innovationlong$trial <- factor(innovationlong$trial, levels = c("Bumpy", "Folded", "Cap1", "Cap2"))

# Descriptive:
hist(innovationlong$time)
# Note that time to solve has a steeply decreasing, long tail distribution -
# suitable for a log transformation
innovationlong$logtime <- log(innovationlong$time)

## Now we have all the data in long format - one bee trial per row - and in 
## wide format - one bee per row, with all the trials.

# We're also going to make a separate 'long' table of abandoning times
abandonedinnovation <- beedata %>%
  pivot_longer(cols = c(F1givinguptime, F2givinguptime, cap1givinguptime, cap2givinguptime), 
               names_to = "trial", 
               values_to = "givinguptime") %>%
  mutate(trial = recode(trial, 
                        "F1givinguptime" = "Bumpy", 
                        "F2givinguptime" = "Folded", 
                        "cap1givinguptime" = "Cap1", 
                        "cap2givinguptime" = "Cap2")) 

## Subsets --------------------------------------------------

# Removing observations where the bee did not visit the flower
innovationlanded <- innovationlong %>%
  filter(landed != 0)

# Make data frame for just the time to solve 
innovationsuccess <- innovationlanded %>%
  filter(solve != 0)


#### SUMMARIES ------------------------------

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

#### EFFECT OF ENVIRONMENT ----------------------------

###### On landing:
glm_landing_fixed <- glm(landed ~ env + trial, family = "binomial", data = innovationlong)
summary(glm_landing_fixed)
# p=0.48 for env
# trial n.s.

###### On solving:
glm_solve_all <- glm(solve ~ env + trial, data = innovationlanded, family = "binomial")
summary(glm_solve_all)
# nothing sign

# FIGURE 5 -----------------------------------

# STEP 1: Summarize landing outcomes
landing_summary <- innovationlong %>%
  filter(!is.na(landed)) %>%
  group_by(env, outcome = ifelse(landed == 1, "Landed", "Did not land")) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(category = "Landing")

# STEP 2: Summarize solving outcomes (only for bees that landed)
solve_summary <- innovationlong %>%
  filter(landed == 1 & !is.na(solve)) %>%
  group_by(env, outcome = ifelse(solve == 1, "Solved", "Did not solve")) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(category = "Solving")

# STEP 3: Combine
plot_data <- bind_rows(landing_summary, solve_summary)

# STEP 4: Set outcome as a factor to control legend order
plot_data$outcome <- factor(
  plot_data$outcome,
  levels = c("Did not land", "Landed", "Did not solve", "Solved")
)

# STEP 5: Define muted custom colors
custom_colors <- c(
  "Did not land" = simplecomplexcolors[1],
  "Landed" = simcomcolors_dark[1],
  "Did not solve" = simplecomplexcolors[2],
  "Solved" = simcomcolors_dark[2]
)

# STEP 6: Plot
sample_sizes <- plot_data %>%
  group_by(env, category) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  mutate(label = paste0("n = ", total))

ggplot(plot_data, aes(x = env, y = count, fill = outcome)) +
  geom_col(position = "fill", width = 0.6) +
  facet_wrap(~ category, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = c("c" = "Complex", "s" = "Simple"))+
  labs(
    x = "Environment",
    y = "Proportion of bees",
    fill = NULL
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.text = element_text(size = 14)
  ) +
  geom_text(
    aes(label = count),
    position = position_fill(vjust = 0.5),
    color = "black",  # <-- black text labels
    size = 4
  )
##!!!! Add asterisk next to 'landing'


##### On time to solve:
env_mod <- lm(logtime ~ env + trial, data = innovationlanded)
summary(env_mod)

# On time to give up:
# Descriptive:
hist(abandonedinnovation$givinguptime)
abandoning_mod <- lm(givinguptime ~ env + trial, 
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

# On landing:
model_SRI_landed <- lm(prop_landed ~ SRI, data = beedata)
model_HB10_landed <- lm(prop_landed ~ HB10, data = beedata)
model_resp_landed <- lm(prop_landed ~ resp, data = beedata)
summary(model_SRI_landed)
summary(model_HB10_landed)
summary(model_resp_landed)
# Again, not a lot of cases of not landing, so not really expecting an effect

# On solving:
model_SRI_solved <- lm(prop_solved ~ SRI, data = beedata)
model_HB10_solved <- lm(prop_solved ~ HB10, data = beedata)
model_resp_solved <- lm(prop_solved ~ resp, data = beedata)
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

#### EFFECT OF FLOWER TYPE ----------------------------
# On time to solve:
trialmod <- aov(logtime ~ trial, data = innovationsuccess)
summary(trialmod)
TukeyHSD(trialmod)

