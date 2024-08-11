#load required libraries
library(reshape2)
library(ggplot2)
library(coxme)
library(survminer)
library(dplyr)
se <- function(x) sd(x)/sqrt(length(x)) #standard error function

#read in data
beeAll <- read.csv(url("https://raw.githubusercontent.com/shannonmcwaters/Innovation/main/BeeInnovationData.csv"))
# Bee's that did failed to visit at least one of the trial flowers
values_to_remove <- c("A10", "A98", "B182", "C32", "Y35", "Y60", "Y48")
#Filter the dataframe to remove rows with matching BeeID values
beeTrimmed <- beeAll %>%
  filter(!BeeID %in% values_to_remove)
################################################################################
#subset data to relevant variables (will add additional variables from main dataset below; for now this simplifies things)
bee3a <- beeAll[,c('BeeID','Env','colony', 'Flower1search','Flower1handling')]
#BeeID: unique id number for each bee
#env: whether the bee was subject to a simple or complex environment
#colony: colony id
#Flower1search: time to alight on first novel flower
#Flower1handling: time spent drinking (head in first novel flower)

bee3a$diff45 <- (beeAll$T5travel - beeAll$T4travel)/beeAll$T4travel #differences between travel times for flowers 4 and 5
bee3a$F1time <- bee3a$Flower1search + bee3a$Flower1handling #total time for 1st novel flower
bee3a$F1solve <- as.numeric(as.factor(beeAll$Flower1Solve))-1 #binary: success on 1st novel flower or not?
bee3a$F2search <- beeAll$Flower2search #time to find 2nd novel flower
bee3a$F2hand <- beeAll$Flower2handling #time to handle 2nd novel flower
bee3a$F2time <- bee3a$F2search + bee3a$F2hand #total time on 2nd novel flower
bee3a$F2solve <- as.numeric(as.factor(beeAll$Flower2Solve))-1 #binary: success on 2nd novel flower
bee3a$cap1hand <- beeAll$Cap1.handling #handling time on 3rd novel flower
bee3a$cap1solve <- as.numeric(as.factor(beeAll$Cap1solve))-1 #binary: success on 3rd novel flower
bee3a$cap2hand <- beeAll$Cap2handling #handling time on 4th novel flower
bee3a$cap2solve <- as.numeric(as.factor(beeAll$Cap2solve))-1 #binary: success on 4th novel flower

#Combined Cox regression considering all flowers together
a <- bee3a[,c('Flower1handling', 'F1solve','Env','diff45')]
names(a) <- c('time','solve','env','resp')
b <- bee3a[,c('F2hand', 'F2solve','Env','diff45')]
names(b) <- c('time','solve','env','resp')
d <- bee3a[,c('cap1hand', 'cap1solve','Env','diff45')]
names(d) <- c('time','solve','env','resp')
e <- bee3a[,c('cap2hand', 'cap2solve','Env','diff45')]
names(e) <- c('time','solve','env','resp')

togeth <- rbind(a,b,d,e)
togeth$BeeID <- rep(1:38,4)
togeth$trial <- rep(letters[22:25], each = 38)

mod <- coxme(Surv(time, solve) ~ env + resp + trial + (1|BeeID), data = togeth)
mod
cox.zph(mod)

mod1 <- coxme(Surv(time, solve) ~ env + resp  + (1|BeeID), data = togeth)
anova(mod,mod1) #test for significant difference between trials

mod2 <- coxph(Surv(time, solve) ~ env + resp  + trial, data = togeth)
anova(mod,mod2) #test for significant individual differences

#Cox regressions for individual trials
summary(m.F1 <- coxph(Surv(F1time,F1solve) ~ Env + diff45, data = bee3a))
cox.zph(m.F1)
summary(m.F2 <- coxph(Surv(F2time,F2solve) ~ Env + diff45, data = bee3a))
cox.zph(m.F2)
summary(m.C1 <- coxph(Surv(cap1hand,cap1solve) ~ Env + diff45, data = bee3a))
cox.zph(m.C1)
summary(m.C2 <- coxph(Surv(cap2hand,cap2solve) ~ Env + diff45, data = bee3a))
cox.zph(m.C2)

#Survival plots
env.dat <- data.frame(Env = c('c','s'),
                      diff45 = rep(mean(bee3a$diff45, na.rm = T), 2))
f1 <- ggsurvplot(survfit(
  coxph(Surv(F1time,F1solve) ~ Env + diff45, data = bee3a),
  data = bee3a, newdata = env.dat),
  conf.int = TRUE,
  legend = 'none',
  legend.labs= c("complex", "simple"),
  palette = c('blue','darkgreen'),
  xlim = c(0,600))

f1 <- f1 + labs(x = "handling  [seconds]", y = "Prop. failing to access sucrose")

f2 <- ggsurvplot(survfit(
  coxph(Surv(F2time,F2solve) ~ Env + diff45, data = bee3a),
  data = bee3a, newdata = env.dat),
  conf.int = TRUE,
  legend = 'none',
  legend.labs=c("complex", "simple"),
  palette = c('blue','darkgreen'),
  xlim = c(0,600))

f2 <- f2 + labs(x = "handling [seconds]", y = "Prop. failing to access sucrose")

c1 <- ggsurvplot(survfit(
  coxph(Surv(cap1hand,cap1solve) ~ Env + diff45, data = bee3a),
  data = bee3a, newdata = env.dat),
  conf.int = TRUE,
  legend = 'none',
  legend.labs=c("complex", "simple"),
  palette = c('blue','darkgreen'),
  xlim = c(0,600))

c1 <- c1 + labs(x = "handling [seconds]", y = "Prop. failing to access sucrose")


bymedian = with(togeth, reorder(BeeID, time, median))
boxplot(time~bymedian,data=togeth,ylim=c(0,60),col=c("grey","black")[as.factor(env)],pch=19)
points(time~bymedian, data=togeth, col=c("blue","green","red","brown")[as.factor(trial)],pch=19)

togeth2 =  subset(togeth, BeeID !=8)
togeth2 %>%
  ggplot() +
  aes(x = time, y = resp, colour = env, group = env) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(trial))+
  geom_smooth(method = "lm")

togeth2 %>%
  ggplot() +
  aes(x = time, y = resp, colour = trial, group = trial) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_hue(direction = 1) +
  theme_minimal()+
  geom_smooth(method = "lm")+
  labs(x="Time to access reward (innovation)", y="Explorative-ness")

ggplot(togeth) +
  aes(x = time, y = resp, colour = env, group = env) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(trial), scales = "free")+
  geom_smooth(method = "lm")+
  labs(x="Time to access reward (innovation)", y="Explorative-ness")

ggplot(togeth) +
  aes(x = time, y = resp, colour = trial, group = trial) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_hue(direction = 1) +
  theme_minimal()+
  geom_smooth(method = "lm")



c2 <- ggsurvplot(survfit(
  coxph(Surv(cap2hand,cap2solve) ~ Env + diff45, data = bee3a),
  data = bee3a, newdata = env.dat),
  conf.int = TRUE,
  legend = 'none',
  legend.labs=c("complex", "simple"),
  palette = c('blue','darkgreen'),
  xlim = c(0,600))
c2 <- c2 + labs(main = "Trial 8b", x = "handling  [seconds]", y = "Prop. failing to access sucrose")

#breaking down the Cox proportional hazards models to see about mechanism
#These models generate Table 1

f1search.g <- glm(Flower1search ~ Env + diff45, family = Gamma(link = 'log'), data = bee3a)
f2search.g <- glm(F2search ~ Env + diff45, family = Gamma(link = 'log'), data = bee3a)

s1 <- bee3a[which(bee3a$Flower1handling > 0),]
f1hand.g <- glm(Flower1handling ~ Env + diff45, family = Gamma(link = 'log'), data = s1)
s2 <- bee3a[which(bee3a$F2hand > 0),]
f2hand.g <- glm(F2hand ~ Env + diff45, family = Gamma(link = 'log'), data = s2)
s3 <- bee3a[which(bee3a$cap1hand > 0),]
cap1hand.g <- glm(cap1hand ~ Env + diff45, family = Gamma(link = 'log'), data = s3)
s4 <- bee3a[which(bee3a$cap2hand > 0),]
cap2hand.g <- glm(cap2hand ~ Env + diff45, family = Gamma(link = 'log'), data = s4)

#Means and standard errors in Table 1
tapply(bee3a$Flower1search, bee3a$Env, length)
tapply(bee3a$Flower1search, bee3a$Env, mean)
tapply(bee3a$Flower1search, bee3a$Env, se)

tapply(bee3a$F2search, bee3a$Env, length)
tapply(bee3a$F2search, bee3a$Env, mean)
tapply(bee3a$F2search, bee3a$Env, se)

tapply(s1$Flower1handling, s1$Env, length)
tapply(s1$Flower1handling, s1$Env, mean)
tapply(s1$Flower1handling, s1$Env, se)

tapply(s2$F2hand, s2$Env, length)
tapply(s2$F2hand, s2$Env, mean)
tapply(s2$F2hand, s2$Env, se)

tapply(s3$cap1hand, s3$Env, length)
tapply(s3$cap1hand, s3$Env, mean)
tapply(s3$cap1hand, s3$Env, se)

tapply(s4$cap2hand, s4$Env, length)
tapply(s4$cap2hand, s4$Env, mean)
tapply(s4$cap2hand, s4$Env, se)
