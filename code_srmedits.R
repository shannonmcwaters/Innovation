#load required libraries
library(reshape2)
library(ggplot2)
library(coxme)
library(survminer)
library(dplyr)
library(emmeans)
se <- function(x) sd(x)/sqrt(length(x)) #standard error function

#read in data
beeAll <- read.csv(url("https://raw.githubusercontent.com/shannonmcwaters/Innovation/main/BeeInnovationData.csv"))
# Bee's that did failed to visit at least one of the trial flowers
values_to_remove <- c("A10", "A98", "B182", "C32", "Y35", "Y60", "Y48","A69")
#Filter the dataframe to remove rows with matching BeeID values
beeTrimmed <- beeAll %>%
  filter(!BeeID %in% values_to_remove)
################################################################################
#subset data to relevant variables (will add additional variables from main dataset below; for now this simplifies things)
bee3 <- beeTrimmed[,c('BeeID','Env','colony', 'Flower1search','Flower1handling')]
#BeeID: unique id number for each bee
#env: whether the bee was subject to a simple or complex environment
#colony: colony id
#Flower1search: time to alight on first novel flower
#Flower1handling: time spent drinking (head in first novel flower)

bee3$diff45 <- (beeTrimmed$T5travel - beeTrimmed$T4travel)/beeTrimmed$T4travel #differences between travel times for flowers 4 and 5
bee3$F1time <- bee3$Flower1handling #total time for 1st novel flower
bee3$F1solve <- as.numeric(as.factor(beeTrimmed$Flower1Solve))-1 #binary: success on 1st novel flower or not?
bee3$F2search <- beeTrimmed$Flower2search #time to find 2nd novel flower
bee3$F2hand <- beeTrimmed$Flower2handling #time to handle 2nd novel flower
bee3$F2time <- bee3$F2hand #total time on 2nd novel flower
bee3$F2solve <- as.numeric(as.factor(beeTrimmed$Flower2Solve))-1 #binary: success on 2nd novel flower
bee3$cap1hand <- beeTrimmed$Cap1.handling #handling time on 3rd novel flower
bee3$cap1solve <- as.numeric(as.factor(beeTrimmed$Cap1solve))-1 #binary: success on 3rd novel flower
bee3$cap2hand <- beeTrimmed$Cap2handling #handling time on 4th novel flower
bee3$cap2solve <- as.numeric(as.factor(beeTrimmed$Cap2solve))-1 #binary: success on 4th novel flower

#Combined Cox regression considering all flowers together
Ta <- bee3[,c('Flower1handling', 'F1solve','Env','diff45')]
names(Ta) <- c('time','solve','env','resp')
Tb <- bee3[,c('F2hand', 'F2solve','Env','diff45')]
names(Tb) <- c('time','solve','env','resp')
Td <- bee3[,c('cap1hand', 'cap1solve','Env','diff45')]
names(Td) <- c('time','solve','env','resp')
Te <- bee3[,c('cap2hand', 'cap2solve','Env','diff45')]
names(Te) <- c('time','solve','env','resp')

Trimtogeth <- rbind(Ta,Tb,Td,Te)
Trimtogeth$BeeID <- rep(1:31,4)
Trimtogeth$trial <- rep(letters[22:25], each = 31)

Trimmedmod <- coxme(Surv(time, solve) ~ env * trial + resp * trial + (1|BeeID), data = Trimtogeth)
Trimmedmod
cox.zph(Trimmedmod)

emmeans_env_by_trial <- emmeans(Trimmedmod, pairwise ~ trial | env, type = "response")

# Test how 'resp' affects each level of 'trial'
emmeans_resp_by_trial <- emmeans(Trimmedmod, pairwise ~ trial | resp, type = "response")

# View the results
emmeans_env_by_trial$contrasts
emmeans_resp_by_trial$contrasts
Trimmedmod1 <- coxme(Surv(time, solve) ~ env + resp  + (1|BeeID), data = Trimtogeth)
anova(Trimmedmod,Trimmedmod1) #test for significant difference between trials

Trimmedmod2 <- coxph(Surv(time, solve) ~ env + resp  + trial, data = Trimtogeth)
anova(Trimmedmod,Trimmedmod2) #test for significant individual differences

#Compare between trials
emmeans_results <- emmeans(Trimmedmod, pairwise ~ trial, type = "response")

#Survival plots
Tenv.dat <- data.frame(Env = c('c','s'),
                      diff45 = rep(mean(bee3$diff45, na.rm = T), 2))
Tf1 <- ggsurvplot(survfit(
  coxph(Surv(F1time,F1solve) ~ Env + diff45, data = bee3),
  data = bee3, newdata = Tenv.dat),
  conf.int = TRUE,
  legend = 'none',
  legend.labs= c("complex", "simple"),
  palette = c('blue','darkgreen'),
  xlim = c(0,600))

Tf1 <- Tf1 + labs(x = "handling  [seconds]", y = "Prop. failing to access sucrose")

Tf2 <- ggsurvplot(survfit(
  coxph(Surv(F2time,F2solve) ~ Env + diff45, data = bee3),
  data = bee3, newdata = Tenv.dat),
  conf.int = TRUE,
  legend = 'none',
  legend.labs=c("complex", "simple"),
  palette = c('blue','darkgreen'),
  xlim = c(0,600))

Tf2 <- Tf2 + labs(x = "handling [seconds]", y = "Prop. failing to access sucrose")

Tc1 <- ggsurvplot(survfit(
  coxph(Surv(cap1hand,cap1solve) ~ Env + diff45, data = bee3),
  data = bee3, newdata = Tenv.dat),
  conf.int = TRUE,
  legend = 'none',
  legend.labs=c("complex", "simple"),
  palette = c('blue','darkgreen'),
  xlim = c(0,600))

Tc1 <- Tc1 + labs(x = "handling [seconds]", y = "Prop. failing to access sucrose")


Tbymedian = with(Trimtogeth, reorder(BeeID, time, median))
boxplot(time~Tbymedian,data=Trimtogeth,ylim=c(0,60),col=c("grey","black")[as.factor(env)],pch=19)
points(time~Tbymedian, data=Trimtogeth, col=c("blue","green","red","brown")[as.factor(trial)],pch=19)

Trimtogeth2 =  subset(Trimtogeth, BeeID !=8)
Trimtogeth2 %>%
  ggplot() +
  aes(x = time, y = resp, colour = env, group = env) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(trial))+
  geom_smooth(method = "lm")

Trimtogeth2 %>%
 ggplot() +
 aes(x = time, y = resp, colour = trial, group = trial) +
 geom_point(shape = "circle", size = 1.5) +
 scale_color_hue(direction = 1) +
 theme_minimal()+
  geom_smooth(method = "lm")+
  labs(x="Time to access reward (innovation)", y="Explorative-ness")

ggplot(Trimtogeth) +
  aes(x = time, y = resp, colour = env, group = env) +
 geom_point(shape = "circle", size = 1.5) +
 scale_color_hue(direction = 1) +
 theme_minimal() +
 facet_wrap(vars(trial), scales = "free")+
  geom_smooth(method = "lm")+
  labs(x="Time to access reward (innovation)", y="Explorative-ness")

ggplot(Trimtogeth) +
 aes(x = time, y = resp, colour = trial, group = trial) +
 geom_point(shape = "circle", size = 1.5) +
 scale_color_hue(direction = 1) +
 theme_minimal()+
  geom_smooth(method = "lm")



Tc2 <- ggsurvplot(survfit(
  coxph(Surv(cap2hand,cap2solve) ~ Env + diff45, data = bee3),
  data = bee3, newdata = Tenv.dat),
  conf.int = TRUE,
  legend = 'none',
  legend.labs=c("complex", "simple"),
  palette = c('blue','darkgreen'),
  xlim = c(0,600))
Tc2 <- Tc2 + labs(main = "Trial 8b", x = "handling  [seconds]", y = "Prop. failing to access sucrose")

#breaking down the Cox proportional hazards models to see about mechanism
#These models generate Table 1

Tf1search.g <- glm(Flower1search ~ Env + diff45, family = Gamma(link = 'log'), data = bee3)
Tf2search.g <- glm(F2search ~ Env + diff45, family = Gamma(link = 'log'), data = bee3)

Ts1 <- bee3[which(bee3$Flower1handling > 0),]
Tf1hand.g <- glm(Flower1handling ~ Env + diff45, family = Gamma(link = 'log'), data = Ts1)
Ts2 <- bee3[which(bee3$F2hand > 0),]
Tf2hand.g <- glm(F2hand ~ Env + diff45, family = Gamma(link = 'log'), data = Ts2)
Ts3 <- bee3[which(bee3$cap1hand > 0),]
Tcap1hand.g <- glm(cap1hand ~ Env + diff45, family = Gamma(link = 'log'), data = Ts3)
Ts4 <- bee3[which(bee3$cap2hand > 0),]
Tcap2hand.g <- glm(cap2hand ~ Env + diff45, family = Gamma(link = 'log'), data = Ts4)

#Means and standard errors in Table 1
tapply(bee3$Flower1search, bee3$Env, length)
tapply(bee3$Flower1search, bee3$Env, mean)
tapply(bee3$Flower1search, bee3$Env, se)

tapply(bee3$F2search, bee3$Env, length)
tapply(bee3$F2search, bee3$Env, mean)
tapply(bee3$F2search, bee3$Env, se)

tapply(Ts1$Flower1handling, Ts1$Env, length)
tapply(Ts1$Flower1handling, Ts1$Env, mean)
tapply(Ts1$Flower1handling, Ts1$Env, se)

tapply(Ts2$F2hand, Ts2$Env, length)
tapply(Ts2$F2hand, Ts2$Env, mean)
tapply(Ts2$F2hand, Ts2$Env, se)

tapply(Ts3$cap1hand, Ts3$Env, length)
tapply(Ts3$cap1hand, Ts3$Env, mean)
tapply(Ts3$cap1hand, Ts3$Env, se)

tapply(Ts4$cap2hand, Ts4$Env, length)
tapply(Ts4$cap2hand, Ts4$Env, mean)
tapply(Ts4$cap2hand, Ts4$Env, se)

