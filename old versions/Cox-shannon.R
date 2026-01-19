rawdata <- read.csv(url("https://raw.githubusercontent.com/shannonmcwaters/Innovation/main/BeeInnovationData.csv"))
####packages####
library(lme4) #needed for GLMMs
library(lmerTest) #needed for obtaining p-values in lmm
library(emmeans) #post-hoc comparisons
library(ggplot2)#for plots
library(dplyr)
library(coxme)
####Data Wrangling####
# Bee's that need to be removed from analysis
bees_to_remove <- c("A69","A10")
#Filter the dataframe to remove rows with matching BeeID values
beeTrimmed <- rawdata %>%
  filter(!BeeID %in% bees_to_remove)
innovatedata <- beeTrimmed[,c('BeeID','Env','colony')]

innovatedata$resp <- (beeTrimmed$T5travel - beeTrimmed$T4travel)/beeTrimmed$T4travel #differences between travel times for flowers 4 and 5
innovatedata$F1time <- beeTrimmed$Flower1handling   #total time for 1st novel flower
innovatedata$F1solve <- as.numeric(as.factor(beeTrimmed$Flower1Solve))-1 #binary: success on 1st novel flower or not?
innovatedata$F2time <- beeTrimmed$Flower2handling #total time on 2nd novel flower
innovatedata$F2solve <- as.numeric(as.factor(beeTrimmed$Flower2Solve))-1 #binary: success on 2nd novel flower
innovatedata$cap1time <- beeTrimmed$Cap1.handling #handling time on 3rd novel flower
innovatedata$cap1solve <- as.numeric(as.factor(beeTrimmed$Cap1solve))-1 #binary: success on 3rd novel flower
innovatedata$cap2time <- beeTrimmed$Cap2handling #handling time on 4th novel flower
innovatedata$cap2solve <- as.numeric(as.factor(beeTrimmed$Cap2solve))-1 #binary: success on 4th novel flower

Bumpy <- innovatedata[,c('F1time', 'F1solve','Env','resp')]
names(Bumpy) <- c('time','solve','env','resp')
Folded <- innovatedata[,c('F2time', 'F2solve','Env','resp')]
names(Folded) <- c('time','solve','env','resp')
Cap1 <- innovatedata[,c('cap1time', 'cap1solve','Env','resp')]
names(Cap1) <- c('time','solve','env','resp')
Cap2 <- innovatedata[,c('cap2time', 'cap2solve','Env','resp')]
names(Cap2) <- c('time','solve','env','resp')

innovatedatalong <- rbind(Bumpy,Folded,Cap1,Cap2)
innovatedatalong$BeeID <- rep(1:36,4)
names = c("bumpy","folded","cap1","cap2")
innovatedatalong$trial <- rep(names, each = 36)

# removing observations where the bee did not visit the flower
innovatedatalong <- innovatedatalong %>%
  filter(time != 0)

# Fit Cox mixed-effects model with random effect for BeeID
coxph(Surv(time, solve) ~ env + trial + BeeID, data = innovatedatalong)
coxme_model <- coxme(Surv(time, solve) ~ env + trial * resp + (1|as.factor(BeeID)), data = innovatedatalong)
coxme_model2 <- coxme(Surv(time, solve) ~ env + trial * resp + (1|as.factor(BeeID)), data = innovatedatalong)

summary(coxme_model)


# Prepare data for hazard ratios
fixed_effects <- data.frame(
  term = c("env", "trialcap1", "trialcap2", "trialfolded", "resp", 
           "trialcap1:resp", "trialcap2:resp", "trialfolded:resp"),
  coef = c(-0.05189, -1.08139, -1.46507, 0.31364, 0.43870, 
           -0.57451, -0.30854, -0.71195)
)

# Calculate hazard ratios and confidence intervals
fixed_effects <- fixed_effects %>%
  mutate(exp_coef = exp(coef),
         lower_ci = exp(coef - 1.96 * sqrt(diag(vcov(coxme_model)))),
         upper_ci = exp(coef + 1.96 * sqrt(diag(vcov(coxme_model)))))

# Hazard ratio plot
ggplot(fixed_effects, aes(x = reorder(term, exp_coef), y = exp_coef)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  coord_flip() +
  labs(x = "Predictors", y = "Hazard Ratio (exp(coef))", 
       title = "Hazard Ratios from Cox Mixed-Effects Model") +
  theme_minimal()

# Interaction plot for responsiveness and trial "cap1"
interaction_data <- innovatedatalong %>%
  group_by(resp, trial) %>%
  summarize(mean_time = mean(time[solve == 1], na.rm = TRUE), .groups = 'drop')

ggplot(innovatedatalong, aes(x = resp, y = time, color = trial)) +
  geom_point(alpha = 0.6, size = 2) +  # Scatter points
  geom_smooth(method = "lm", se = F, aes(group = trial), size = 1) +  # Linear trend lines
  labs(title = "Time to Solve vs. Responsiveness by Trial Type",
       x = "Responsiveness",
       y = "Time to Solve",
       color = "Trial Type") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set1")  # Color palette for clarity
