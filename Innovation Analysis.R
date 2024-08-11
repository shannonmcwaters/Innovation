rawdat <- read.csv(url("https://raw.githubusercontent.com/shannonmcwaters/Innovation/main/BeeInnovationData.csv"))
####packages####
library(lme4) #needed for GLMMs
library(lmerTest) #needed for obtaining p-values in lmm
library(emmeans) #post-hoc comparisons
library(ggplot2)#for plots
####Data Wrangling####
# Bee's that need to be removed from analysis
bees_to_remove <- c("A69")
#Filter the dataframe to remove rows with matching BeeID values
beeTrimmed <- rawdat %>%
  filter(!BeeID %in% bees_to_remove)

innovatedata <- beeTrimmed[,c('BeeID','Env','colony')]
#BeeID: unique id number for each bee
#env: whether the bee was subject to a simple or complex environment
#colony: colony id
#Flower1search: time to alight on first novel flower
#Flower1handling: time spent drinking (head in first novel flower)

innovatedata$resp <- (beeTrimmed$T5travel - beeTrimmed$T4travel)/beeTrimmed$T4travel #differences between travel times for flowers 4 and 5
innovatedata$F1time <- beeTrimmed$Flower1handling #total time for 1st novel flower
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
innovatedatalong$BeeID <- rep(1:37,4)
innovatedatalong$trial <- rep(letters[22:25], each = 37)
# removing observations where the bee did not visit the flower
innovatedatalong <- innovatedatalong %>%
  filter(time != 0)
#data of all bees
solve = innovatedatalong
#make data frame for just the time to solve 
innovation <- innovatedatalong %>%
  filter(solve != 0)
#make data frame for just the time to give up
abandoning <- innovatedatalong %>%
  filter(solve != 1)


####Did the bee solve?####
solve_mod<- glmer(solve ~ env + resp * trial + (1 | BeeID), 
                    data = solve, 
                    family = binomial)
summary(solve_mod)
# Obtain estimated marginal means (EMMs) for 'trial'
solvemeans <- emmeans(solve_mod, ~ trial)
# Perform pairwise comparisons between trials
pairwise_comparisons <- contrast(solvemeans, method = "pairwise")
# Summary of pairwise comparisons
summary(pairwise_comparisons)

##plots!##
# Calculate the proportion of bees that solved for each trial
proportion_summary <- solve %>%
  group_by(trial) %>%
  summarize(proportion_solved = mean(solve == 1), .groups = 'drop')

# Plot the proportion of solved cases for each trial
ggplot(heatmap_data, aes(x = trial, y = BeeID, fill = factor(solved))) +
  geom_tile(color = "white") +  # Add white lines between tiles
  scale_fill_manual(values = c("0" = "black", "1" = "grey"), labels = c("Not Solved", "Solved")) +  # Black and grey colors
  scale_x_discrete(labels = c("bumpy", "folded", "cap1", "cap2")) +  # Ensure x-axis labels
  labs(x = "Trial", y = "BeeID", fill = "Solved") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_blank())  # Remove plot title
####Time to solve####
innovation$logtime <- log(innovation$time)
innovation_mod <- lmer(logtime ~ env + trial * resp + (1 | BeeID), 
                       data = innovation)
summary(innovation_mod)
isSingular(innovation_mod)
innovation_mod2 <- lm(logtime ~ env + resp * trial,
                       data = innovation)
summary(innovation_mod2)

# Obtain estimated marginal means (EMMs) for 'trial'
innovmeans <- emmeans(innovation_mod, ~ trial)
# Perform pairwise comparisons between trials
pairwise_comparisons2 <- contrast(innovmeans, method = "pairwise")
# Summary of pairwise comparisons
summary(pairwise_comparisons2)

#Plot ave times per trial#
ggplot(innovation, aes(x = trial, y = time, fill = env)) +
  geom_boxplot() +
  scale_fill_manual(values = c("c" = "black", "s" = "grey"), 
                    labels = c("complex", "simple")) +  # Update legend labels
  labs(x = "Trial", y = "Time to solve (seconds)", fill = "Environment") +  # Change the legend title to "Type"
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = c("bumpy", "folded", "cap1", "cap2"))

#look at trials separately with responsiveness#
average_time_data <- innovation %>%
  group_by(trial, resp, env) %>%
  summarize(avg_time = mean(time, na.rm = TRUE), .groups = 'drop')

# Plot with a 2x2 matrix of scatterplots and lines of best fit
ggplot(average_time_data, aes(x = resp, y = avg_time, color = env)) +
  geom_point(size = 3) +  # Add points for average time
  geom_smooth(method = "lm", se = FALSE, aes(color = env)) +  # Add line of best fit without confidence interval
  labs(x = "Responsiveness", y = "Time to solve (seconds)", color = "Environment") +
  facet_wrap(~ trial, ncol = 2, labeller = as_labeller(c("v" = "bumpy", "w" = "folded", "x" = "cap1", "y" = "cap2"))) +  # Create a 2x2 matrix of plots, with 2 columns
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###Time to give up####
hist(abandoning$time)
abandoning_mod <- lm(time ~ env + resp + trial, 
                       data = abandoning)
summary(abandoning_mod)
#removing trial as the trials are unevenly represented and there is only one observation for trial v
abandoning_mod2 <- lm(time ~ env + resp,
                       data = abandoning)
summary(abandoning_mod2)

#plot
ggplot(abandoning, aes(x = trial, y = time, color = env)) +
  geom_point(aes(shape = env), size = 3) +  # Add points with different shapes for each env
  labs(x = "Trial", y = "Time", color = "Environment", shape = "Environment") +
  scale_color_manual(values = c("c" = "black", "s" = "grey"), labels = c("complex", "simple")) +  # Update color legend
  scale_shape_manual(values = c("c" = 16, "s" = 17), labels = c("complex", "simple")) +  # Update shape legend
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("bumpy", "folded", "cap1", "cap2"))
