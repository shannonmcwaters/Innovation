# Shannon McWaters, James J. Kearsley, David W. Kikuchi, 
# Timothy J. Polnaszek, Anna Dornhaus
# R script for data analysis and figures for 
# Innovation project: bumble bees, inherent traits, and environment

# Paper reference: XXX (to be updated on acceptance of the paper)
### PAPER VERSION    DRAFT

### Libraries ------------------
library(lme4) #needed for GLMMs
library(lmerTest) #needed for obtaining p-values in lmm
library(emmeans) #post-hoc comparisons
library(coxme)
library(stringdist)
library(tidyverse)
library(scales)  # for number_format
library(sjPlot) # for nice model output tables
library(ggh4x) #for formatting individual x-axes in grid plot
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
beedata$F1search <- beedata$Flower1search
beedata$F2search <- beedata$Flower2search
beedata$tot_search12 <- beedata$F1search + beedata$F2search
# Calculate average time spent traveling over trips 2-4 in the pre-training 
beedata$avgtravel <- rowMeans(beedata[,c("T2_T","T3_T","T4_T")], na.rm=TRUE)
#  resp: responsiveness
beedata$resp <- (beedata$T5travel - beedata$avgtravel)/beedata$avgtravel #differences between travel times for trip 5 and avg of trips before
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

# Removing columns we don't need or that repeat information
# Long format: one row per trial (=trip), multiple per bee
innovationlong <- innovationlong %>% 
  select(
    BeeID,
    colony,
    env,
    trial, 
    landed, solve,
    time,
    logtime,
    SRI, resp, HB10, H_F1_T1, tot_search12, avgtravel
  )
abandonedinnovation <- abandonedinnovation %>% 
  select(
    BeeID,
    colony,
    env,
    trial, 
    givinguptime,
    SRI, resp, HB10, H_F1_T1, tot_search12, avgtravel
  )

# Wide format: one row per bee, with columns for each trait
# We don't need information for each trial any more for this dataset.
beedata <- beedata %>%
  select(
    BeeID,
    colony,
    env,
    prop_landed, prop_solved, 
    avg_time,
    SRI, resp, HB10, H_F1_T1, tot_search12, avgtravel
  )

# Removing observations where the bee did not visit the flower
innovationlanded <- innovationlong %>%
  filter(landed != 0)

# Make data frame for just the time to solve 
innovationsuccess <- innovationlanded %>%
  filter(solve != 0)

# Descriptive (Fig S1):
ggplot(innovationsuccess
       , aes(x = time)
) + 
  geom_histogram(
    breaks=seq(0.1, 60.1, by = 5),
    fill = "grey80",
    color = "black"
  ) +
  geom_vline(
    xintercept = median(innovationlong$time[innovationlong$time > 0], na.rm = TRUE),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_x_continuous(
    breaks = seq(0, 60, by = 10)
  ) +
  labs(
    x = "Time to solve [s]",
    y = "Frequency"
  ) +
  theme_classic(base_size = 14)
# Note that time to solve has a steeply decreasing, long tail distribution -
# suitable for a log transformation


#### SUMMARIES ------------------------------

# Calculate mean solving time for each trial type
mean_times <- innovationsuccess %>%
  group_by(trial) %>%
  summarise(mean_time = mean(time, na.rm = TRUE))

# Another version of per-bee (wide) data table, with traits
# Reshape per-bee values for faceted plot
bee_longplot <- beedata %>%
  pivot_longer(cols = c(H_F1_T1, tot_search12, avgtravel, SRI, resp, HB10), names_to = "trait", values_to = "score") %>%
  pivot_longer(cols = c(prop_landed, prop_solved), names_to = "outcome", values_to = "prop") %>%
  mutate(trait = recode(trait,
                        H_F1_T1 = "First handling",
                        tot_search12 = "Search time",
                        avgtravel = "Travel time",
                        SRI = "Routine",
                        resp = "Responsiveness",
                        HB10 = "Exploration"
                        )
         , outcome = recode(outcome,
                            prop_landed = "Prop. trials landed",
                            prop_solved = "Prop. trials solved")
  ) %>%
  mutate(trait = factor(trait, levels = c("First handling",
                                          "Search time", 
                                          "Travel time",
                                          "Routine",
                                          "Responsiveness",
                                          "Exploration"))
         , outcome = factor(outcome, levels = c("Prop. trials landed", "Prop. trials solved"))
  )

meanvals <- c(mean(beedata$H_F1_T1)
                , mean(beedata$tot_search12)
                , mean(beedata$avgtravel)
                , mean(beedata$SRI)
                , mean(beedata$resp)
                , mean(beedata$HB10)
)
sdvals <- c(sd(beedata$H_F1_T1)
              , sd(beedata$tot_search12)
              , sd(beedata$avgtravel)
              , sd(beedata$SRI)
              , sd(beedata$resp)
              , sd(beedata$HB10)
)
#bee_longplot$difftomean <- (bee_longplot$score 
#                            - meanvals[bee_longplot$trait])
bee_longplot$zscore <- (bee_longplot$score - meanvals[bee_longplot$trait])/sdvals[bee_longplot$trait]

## ANALYSES -------------------------------------------

#### EFFECT OF ENVIRONMENT ----------------------------

###### On landing: Table S1 ------------------------------------
glm_landing <- glm(landed ~ env + trial, family = "binomial", data = innovationlong)
summary(glm_landing)
# p=0.48 for env
# trial n.s.

# Output table for supplementary
tab_model(glm_landing
          , show.re.var = TRUE
          , pred.labels = c("Intercept",
                          "Environment (complex vs simple)",
                          "Trial (Folded vs Bumpy)",
                          "Trial (Cap1 vs Bumpy)",
                          "Trial (Cap2 vs Bumpy)"
                          )
          , dv.labels = "Effect on landing on novel flower"
)

###### On solving: Table S2 ------------------------------------
glm_solving <- glm(solve ~ env + trial, data = innovationlanded, family = "binomial")
summary(glm_solving)
# nothing sign
tab_model(glm_solving
          , show.re.var = TRUE
          , pred.labels = c("Intercept",
                            "Environment (complex vs simple)",
                            "Trial (Folded vs Bumpy)",
                            "Trial (Cap1 vs Bumpy)",
                            "Trial (Cap2 vs Bumpy)"
          )
          , dv.labels = "Effect on solving novel flower"
)



# FIGURE 3 -----------------------------------

# STEP 1: Summarize landing outcomes
landing_summary <- innovationlong %>%
  filter(!is.na(landed)) %>%
  group_by(env, outcome = ifelse(landed == 1, "Landed", "Did not land")) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(category = "Landing  *")

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
    y = "Proportion of observations",
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


# On time to solve: Table S3 -------------------------------
lm_time_to_solve <- lm(logtime ~ env * trial, data = innovationsuccess)
summary(lm_time_to_solve)
tab_model(lm_time_to_solve
          , show.re.var = TRUE
          , pred.labels = c("Intercept",
                            "Environment (complex vs simple)",
                            "Trial (Folded vs Bumpy)",
                            "Trial (Cap1 vs Bumpy)",
                            "Trial (Cap2 vs Bumpy)",
                            "Env x Trial (Folded)",
                            "Env x Trial (Cap1)",
                            "Env x Trial (Cap2)"
          )
          , dv.labels = "Effect on time to solve (reach reward) on novel flower"
)

# FIGURE 4 --------------------------------------------

# Boxplot of innovation time for each trial
graph_data <- innovationsuccess
ymax <- max(graph_data$time)
offset <- 0.2
N_s <- table(subset(graph_data, env=="s")$trial)
N_c <- table(subset(graph_data, env=="c")$trial)

ggplot(graph_data, aes(x = trial, y = time, fill = env)) +
  geom_boxplot(coef = Inf, position = position_dodge(width = 0.75)) +
  scale_y_log10() +
  scale_fill_manual(
    values = c("s" = simplecomplexcolors[1], "c" = simplecomplexcolors[2]),
    labels = c("Simple", "Complex"),
    name = "Environment"
  )+
  theme_minimal() +
  labs(x = "Trial", y = "Time to Solve [s]") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) +
  annotate("text", x = 1, y = ymax + 2, label = "A", size = 5, color = "black") +
  annotate("text", x = 2, y = ymax + 2, label = "A", size = 5, color = "black") +
  annotate("text", x = 3, y = ymax + 2, label = "B", size = 5, color = "black") +
  annotate("text", x = 4, y = ymax + 2, label = "B", size = 5, color = "black") +
  annotate("text", x = 1+offset, y = 0.9, label = N_c[[1]], size = 5, color = simcomcolors_dark[2]) +
  annotate("text", x = 2+offset, y = 0.9, label = N_c[[2]], size = 5, color = simcomcolors_dark[2]) +
  annotate("text", x = 3+offset, y = 0.9, label = N_c[[3]], size = 5, color = simcomcolors_dark[2]) +
  annotate("text", x = 4+offset, y = 0.9, label = N_c[[4]], size = 5, color = simcomcolors_dark[2]) +
  annotate("text", x = 1-offset, y = 0.9, label = N_s[[1]], size = 5, color = simcomcolors_dark[1]) +
  annotate("text", x = 2-offset, y = 0.9, label = N_s[[2]], size = 5, color = simcomcolors_dark[1]) +
  annotate("text", x = 3-offset, y = 0.9, label = N_s[[3]], size = 5, color = simcomcolors_dark[1]) +
  annotate("text", x = 4-offset, y = 0.9, label = N_s[[4]], size = 5, color = simcomcolors_dark[1]) 


# On time to give up: Table S4 ---------------------------------------------
lm_abandoning <- lm(givinguptime ~ env + trial, 
                       data = abandonedinnovation)
summary(lm_abandoning)
tab_model(lm_abandoning
          , show.re.var = TRUE
          , pred.labels = c("Intercept",
                            "Environment (complex vs simple)",
                            "Trial (Folded vs Bumpy)",
                            "Trial (Cap1 vs Bumpy)",
                            "Trial (Cap2 vs Bumpy)"
          )
          , dv.labels = "Effect on the time to abandon a flower without solving"
)

# FIGURE S2 --------------------------------------------------
graph_data <- subset(abandonedinnovation, givinguptime >0)
graph_data$trial_factor <- factor(graph_data$trial, levels = c("Bumpy", "Folded", "Cap1", "Cap2"))
N_s <- table(subset(graph_data, env=="s")$trial_factor)
N_c <- table(subset(graph_data, env=="c")$trial_factor)

ggplot(graph_data, aes(x = trial_factor,
                       y = givinguptime, fill = factor(env, levels = c("s", "c")))) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  scale_fill_manual(
    values = c("s" = simplecomplexcolors[1], "c" = simplecomplexcolors[2]),
    labels = c("Simple", "Complex"),
    name = "Environment"
  ) +
  theme_minimal() +
  labs(x = "Trial", y = "Time to give up [s]") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) +
  annotate("text", x = 1+offset, y = - 2, label = N_c[[1]], size = 5, color = simcomcolors_dark[2]) +
  annotate("text", x = 2+offset, y = - 2, label = N_c[[2]], size = 5, color = simcomcolors_dark[2]) +
  annotate("text", x = 3+offset, y = - 2, label = N_c[[3]], size = 5, color = simcomcolors_dark[2]) +
  annotate("text", x = 4+offset, y = - 2, label = N_c[[4]], size = 5, color = simcomcolors_dark[2]) +
  annotate("text", x = 1-offset, y = - 2, label = N_s[[1]], size = 5, color = simcomcolors_dark[1]) +
  annotate("text", x = 2-offset, y = - 2, label = N_s[[2]], size = 5, color = simcomcolors_dark[1]) +
  annotate("text", x = 3-offset, y = - 2, label = N_s[[3]], size = 5, color = simcomcolors_dark[1]) +
  annotate("text", x = 4-offset, y = - 2, label = N_s[[4]], size = 5, color = simcomcolors_dark[1]) 


# On other bee traits: -------------------------------
wilcox.test(H_F1_T1 ~ env, data = beedata)
wilcox.test(tot_search12 ~ env, data = beedata)
wilcox.test(avgtravel ~ env, data = beedata)
wilcox.test(SRI ~ env, data = beedata)
wilcox.test(resp ~ env, data = beedata)
wilcox.test(HB10 ~ env, data = beedata)

# FIGURE 5 -----------------------------------
graph_data <- bee_longplot
ggplot(graph_data, aes(x = trait, y = zscore, fill = env)) +
  geom_boxplot(coef = Inf
               , position = position_dodge(width = 0.75)
               , alpha = 0.7
               ) +
  labs(
    x = "Environment",
    y = "Relative trait score [(x-mean)/stdev]"
  ) +
  scale_fill_manual(
    values = simplecomplexcolors,
    labels = c("s" = "Simple", "c" = "Complex"),
    name = "Environment"
  ) +
  scale_x_discrete(labels = c("s" = "Simple", "c" = "Complex")) +
  theme_minimal(base_size = 16) +
  theme(
    strip.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  ) 

#### EFFECT OF FLOWER TYPE ----------------------------
# On time to solve:
trialmod <- aov(logtime ~ trial, data = innovationsuccess)
summary(trialmod)
TukeyHSD(trialmod)

#### EFFECT OF BEE TRAITS ----------------------------

# For landing (yes/no) and solving (yes/no) ---------------
# We use one data point per bee, the proportion of time landing/solving. 

# On landing:
model_firsthandl_landed <- lm(prop_landed ~ H_F1_T1, data = beedata)
model_searchtime_landed <- lm(prop_landed ~ tot_search12, data = beedata)
model_travel <- lm(prop_landed ~ avgtravel, data = beedata)
model_SRI_landed <- lm(prop_landed ~ SRI, data = beedata)
model_resp_landed <- lm(prop_landed ~ resp, data = beedata)
model_HB10_landed <- lm(prop_landed ~ HB10, data = beedata)

summary(model_firsthandl_landed)
summary(model_travel)
summary(model_searchtime_landed) # only one sign
summary(model_SRI_landed)
summary(model_resp_landed)
summary(model_HB10_landed)

# On solving:
model_firsthandl_solved <- lm(prop_solved ~ H_F1_T1, data = beedata)
model_searchtime_solved <- lm(prop_solved ~ tot_search12, data = beedata)
model_travel_solved <- lm(prop_solved ~ avgtravel, data = beedata)
model_SRI_solved <- lm(prop_solved ~ SRI, data = beedata)
model_resp_solved <- lm(prop_solved ~ resp, data = beedata)
model_HB10_solved <- lm(prop_solved ~ HB10, data = beedata)

summary(model_firsthandl_solved)
summary(model_searchtime_solved) # only one sign
summary(model_travel_solved)
summary(model_SRI_solved)
summary(model_resp_solved)
summary(model_HB10_solved)

# FIGURE S3 ------------------------------------------
ggplot(bee_longplot, aes(x = score, y = prop, color = env)) +
  geom_point(alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = alpha("black", 0.2)) +
  facet_grid(
    rows = vars(outcome),
    cols = vars(trait),
    scales = "free_x",
    switch = "y"
  ) +
  ggh4x::facetted_pos_scales(
    x = list(
      # FORCE 0 to appear + axis starts at 0 for this facet
      trait == "First handling" ~ scale_x_continuous(
        limits = c(0, NA),
        breaks = c(0, 5, 10, 15, 20),
        labels = label_number(accuracy = 1)
      ),
      
      trait == "Search Time" ~ scale_x_continuous(
        breaks = pretty_breaks(n = 3),
        labels = label_number(big.mark = ",", accuracy = 1)
      ),
      
      trait == "Routine" ~ scale_x_continuous(
        breaks = c(0.3, 0.5, 0.7),
        labels = label_number(accuracy = 0.1)
      ),
      
      trait == "Responsiveness" ~ scale_x_continuous(
        breaks = c(-0.5, 0, 0.5, 1.0),
        labels = label_number(accuracy = 0.1)
      ),
      
      trait == "Exploration" ~ scale_x_continuous(
        breaks = c(0, 0.5, 1),
        labels = label_number(accuracy = 0.1)
      )
    )
  ) +
  labs(x = "Trait score", y = NULL) +
  theme_minimal(base_size = 14) +
  scale_color_manual(
    values = c("s" = simcomcolors_dark[1], "c" = simcomcolors_dark[2]),
    labels = c("s" = "Simple", "c" = "Complex"),
    name = "Environment"
  ) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, size = 14),
    strip.text.x = element_text(size = 10),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


# On time to solve: --------------------------------------------
# For the continuous outcome of solving time, we use all innovation trials 
# (flower types) separately, thus multiple measurements per bee; 
# but only trials in which the bee actually solved (i.e. accessed reward).

# Handling time of first flower: Table S5 ------------------------
hand_lm <- lm(logtime ~ log(H_F1_T1) * trial, data = innovationsuccess)
summary(hand_lm)
# Close to sign
tab_model(hand_lm
          , show.re.var = TRUE
          , pred.labels = c("Intercept",
                            "First handling time",
                            "Trial (Folded vs Bumpy)",
                            "Trial (Cap1 vs Bumpy)",
                            "Trial (Cap2 vs Bumpy)",
                            "First hand x Trial (Folded)",
                            "First hand x Trial (Cap1)",
                            "First hand x Trial (Cap2)"
          )
          , dv.labels = "Effect on solving time"
)

# Search time in first two innovation trials: Table S6 ---------------
searchmod <- lm(logtime ~ log(tot_search12) * trial, data = innovationsuccess)
summary(searchmod)
tab_model(searchmod
          , show.re.var = TRUE
          , pred.labels = c("Intercept",
                            "Search time on B&F",
                            "Trial (Folded vs Bumpy)",
                            "Trial (Cap1 vs Bumpy)",
                            "Trial (Cap2 vs Bumpy)",
                            "Search x Trial (Folded)",
                            "Search x Trial (Cap1)",
                            "Search x Trial (Cap2)"
          )
          , dv.labels = "Effect on solving time"
)

# Travel time in pretrials: Table S7 ---------------
travelmod <- lm(logtime ~ log(avgtravel) * trial, data = innovationsuccess)
summary(travelmod)
tab_model(travelmod
          , show.re.var = TRUE
          , pred.labels = c("Intercept",
                            "Travel in trial 2-4",
                            "Trial (Folded vs Bumpy)",
                            "Trial (Cap1 vs Bumpy)",
                            "Trial (Cap2 vs Bumpy)",
                            "Travel x Trial (Folded)",
                            "Travel x Trial (Cap1)",
                            "Travel x Trial (Cap2)"
          )
          , dv.labels = "Effect on solving time"
)

# SRI (routine formation): Table S8 ------------------------------
sri_lm <- lm(logtime ~ SRI * trial, data = innovationsuccess)
summary(sri_lm)
# Sign SRI and interaction SRI x trialCap1
tab_model(sri_lm
          , show.re.var = TRUE
          , pred.labels = c("Intercept",
                            "SRI (Routine formation)",
                            "Trial (Folded vs Bumpy)",
                            "Trial (Cap1 vs Bumpy)",
                            "Trial (Cap2 vs Bumpy)",
                            "SRI x Trial (Folded)",
                            "SRI x Trial (Cap1)",
                            "SRI x Trial (Cap2)"
          )
          , dv.labels = "Effect on solving time"
)

# Responsiveness: Table S9 -----------------------------
resp_lm <- lm(logtime ~ resp * trial, data = innovationsuccess)
summary(resp_lm)
# Sign interaction resp x trialFolded
tab_model(resp_lm
          , show.re.var = TRUE
          , pred.labels = c("Intercept",
                            "Responsiveness",
                            "Trial (Folded vs Bumpy)",
                            "Trial (Cap1 vs Bumpy)",
                            "Trial (Cap2 vs Bumpy)",
                            "resp x Trial (Folded)",
                            "resp x Trial (Cap1)",
                            "resp x Trial (Cap2)"
          )
          , dv.labels = "Effect on solving time"
)

# Exploration: Table S10 -----------------------------
exp_lm <- lm(logtime ~ HB10 * trial, data = innovationsuccess)
summary(exp_lm)
# Not sign
tab_model(exp_lm
          , show.re.var = TRUE
          , pred.labels = c("Intercept",
                            "Exploration",
                            "Trial (Folded vs Bumpy)",
                            "Trial (Cap1 vs Bumpy)",
                            "Trial (Cap2 vs Bumpy)",
                            "Expl x Trial (Folded)",
                            "Expl x Trial (Cap1)",
                            "Expl x Trial (Cap2)"
          )
          , dv.labels = "Effect on solving time"
)


# FIGURE 6 ---------------------------------------------
ggplot(bee_longplot, aes(x = score, y = avg_time
#                         , color = env
                         )) +
  geom_smooth(method = "lm", se = TRUE
              , color = alpha("black", 0.2)
              , fill = alpha("grey", 0.05)) +
  geom_point(color = "darkblue", alpha = 0.2) +
  facet_wrap(~ trait, scales = "free_x", nrow = 1) +
  scale_y_log10() +
  scale_color_manual(
#    values = c("s" = simcomcolors_dark[1], 
#               "c" = simcomcolors_dark[2]),
#    labels = c("Simple", "Complex"),
#    name = "Environment"
  ) +
  labs(
    x = "Trait score",
    y = "Mean solving time [s]"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title.y = element_text(size = 22),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


# Bee IDs ----------------------------------
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
# Repeatability is zero



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
# Cap1 & Cap2 are different but env is not sign

