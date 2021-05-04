#Project 1
library(ggplot2)
library(GGally)
library(ggpubr)

plasma <- read.delim("Data/plasma.txt")
head(plasma)
summary(plasma)
# At least one data point was 0, this needs to be taken 
PositivePlasma <- plasma[plasma$betaplasma > 0, ]
summary(PositivePlasma)
#Turn categorical variables into factor variables:
PositivePlasma$sex <- factor(PositivePlasma$sex, 
                               levels = c(1, 2),
                               labels = c("Male", "Female"))

PositivePlasma$smokstat <- factor(PositivePlasma$smokstat,
                                    levels = c(1, 2, 3),
                                    labels = c("Never", "Former", "Current"))

PositivePlasma$bmicat <- factor(PositivePlasma$bmicat,
                                  levels = c(1,2,3,4),
                                  labels = c("Underweight", "Normal", "Overweight", "Obese"))

PositivePlasma$vituse <- factor(PositivePlasma$vituse,
                                  levels = c(1, 2, 3),
                                  labels = c("Yes, fairly often", "Yes, not often", "No"))

#relevel to suitable reference categories:
PositivePlasma$sex<-relevel(PositivePlasma$sex, "Female")
PositivePlasma$smokstat<-relevel(PositivePlasma$smokstat, "Never")
PositivePlasma$bmicat <- relevel(PositivePlasma$bmicat,"Normal")
PositivePlasma$vituse <- relevel(PositivePlasma$vituse,"Yes, fairly often")

#Part 0 a)
a <- 225.498
#Part 0 b) beta-carotene is considered low if below 225.498 ng/ml
PositivePlasma$lowplasma <- as.numeric(PositivePlasma$betaplasma < a)
# factor version for ease of use later on:
PositivePlasma$plasmacat <- factor(PositivePlasma$lowplasma,
                           levels = c(0, 1),
                           labels = c("high", "low"))
# How many are low?
low_total <- sum(PositivePlasma$lowplasma)
high_total <- 315 - low_total
low_percentage <- (low_total/316)*100 

table(PositivePlasma$smokstat, PositivePlasma$lowplasma)


(model.1 <- glm(lowplasma ~ smokstat, family = "binomial", data = PositivePlasma))
summary(model.1)

# beta: log-odds(ratio) with c.i.:
model.1$coefficients
(ci.beta <- confint(model.1))

# Odds (exp(beta0)) and OR, odds ratio, exp(beta1), exp(beta2)
exp(model.1$coefficients)
(ci.or <- exp(ci.beta))

#probability for "Never" to have low beta-plasma:
(prob_0 <- 2.250000/(2.250000+1))
#confidence interval:
(prob_0_low <- 1.6120554/(1.6120554+1))
(prob_0_high <- 3.187487/(3.187487+1))

#probability for "Former" to have low beta-plasma:
(prob_1 <- 1.318008*prob_0/(1-prob_0)/(1+1.318008*prob_0/(1-prob_0)))
#confidence interval:
(prob_1_low <- 0.7707618*prob_0/(1-prob_0)/(1+0.7707618*prob_0/(1-prob_0)))
(prob_1_high <- 2.280987*prob_0/(1-prob_0)/(1+2.280987*prob_0/(1-prob_0)))

#probability for "Current" to have low beta-plasma:
(prob_1 <- 5.925926*prob_0/(1-prob_0)/(1+5.925926*prob_0/(1-prob_0)))
#confidence interval:
(prob_2_low <- 2.0209273*prob_0/(1-prob_0)/(1+2.0209273*prob_0/(1-prob_0)))
(prob_2_high <- 25.326220*prob_0/(1-prob_0)/(1+25.326220*prob_0/(1-prob_0)))

#Wald test
summary(model.1)$coefficients
# Since |0.2761213| < lambda_0.025 = 1.96 we accept
# H0: beta_1 = 0
# Since |1.7793369| < lambda_0.025 = 1.96 we accept
# H0: beta_2 = 0
#The smoking has no significant impact on the beta-carotene levels

# 1b)
(age.model <- glm(lowplasma ~ age, family = "binomial", data = PositivePlasma))

PositivePlasma.pred <- cbind(
  PositivePlasma,
  phat = predict(age.model, type = "response"),
  logit = predict(age.model, se.fit = TRUE))

(Age.plot<-ggplot(PositivePlasma.pred, aes(age, lowplasma)) +
    geom_point() +
    geom_smooth(se = FALSE, linetype = "dashed") +
    geom_line(aes(y = phat), color = "red", size = 1) +
    xlab("age") +
    ylab("probability of having low beta-carotene") +
    labs(caption = "red = fitted line, blue dashed = moving average") +
    theme(text = element_text(size = 14)))

(BetaA.lm <- cbind(summary(age.model)$coefficients,ci =confint(age.model)))

# beta: log-odds(ratio) with c.i.:
age.model$coefficients
(ci.beta <- confint(age.model))

# Odds (exp(beta0)) and OR, odds ratio, exp(beta1)
exp(age.model$coefficients)
(ci.or <- exp(ci.beta))
#Wald test
summary(age.model)$coefficients
# Since |-0.02835726| < lambda_0.025 = 1.96 we can reject
# H0: beta_1 = 0
#The age has a significant impact on the beta-carotene levels

#The odds of having low beta plasma as you age by 1 year decrease by 2.8%
((1-exp(age.model$coefficients[2]))*100)

#The odds of having low beta plasma as you age by 10 years decrease by 24.7% with a 37% decrease lwr and a 10% decrease upr.
((1-exp(10*BetaA.lm[2,c(1,5,6)])))
#Predicted probabilities
(lambda <- qnorm(1 - 0.05/2))
PositivePlasma.pred$logit.lwr <- PositivePlasma.pred$logit.fit - lambda*PositivePlasma.pred$logit.se.fit
PositivePlasma.pred$logit.upr <- PositivePlasma.pred$logit.fit + lambda*PositivePlasma.pred$logit.se.fit
# transform the log-odds intervals into C.I. for odds####
PositivePlasma.pred$odds.lwr <- exp(PositivePlasma.pred$logit.lwr)
PositivePlasma.pred$odds.upr <- exp(PositivePlasma.pred$logit.upr)
# transform the odds intervals into C.I. for p####
PositivePlasma.pred$p.lwr <- PositivePlasma.pred$odds.lwr/(1 + PositivePlasma.pred$odds.lwr)
PositivePlasma.pred$p.upr <- PositivePlasma.pred$odds.upr/(1 + PositivePlasma.pred$odds.upr)

(Q1b.plot<-Age.plot+geom_ribbon(data=PositivePlasma.pred, aes(ymin = p.lwr, ymax = p.upr), alpha = 0.2)+labs(caption = "Red = fitted line, with 95% confidence interval. Blue dashed = moving average"))
#30 & 31 year old
#Prediction intervals based on age
x0<-data.frame(age = 30)
logit0<-round(predict(age.model, newdata = x0, interval = "prediction"),digits = 12)
predint0<-cbind(x0,exp(logit0)/(1+exp(logit0)))

x1<-data.frame(age = 70)
logit1<-round(predict(age.model, newdata = x1, interval = "prediction"),digits = 12)
predint1<-cbind(x1,exp(logit1)/(1+exp(logit1)))
#31 year old & 71 year old
x2<-data.frame(age = 31)
logit2<-round(predict(age.model, newdata = x2, interval = "prediction"),digits = 12)
predint2<-cbind(x2,exp(logit2)/(1+exp(logit2)))

x3<-data.frame(age = 71)
logit3<-round(predict(age.model, newdata = x3, interval = "prediction"),digits = 12)
predint3<-cbind(x3,exp(logit3)/(1+exp(logit3)))


#Differences: what is the slope and what are the reasons for the differences between the 2 age groups
(predint2[2]-predint0[2])#Difference btwn 30 & 31: 0.00377553
(predint3[2]-predint1[2])#Difference btwn 70 & 71: 0.006600449

#Question c:
PositivePlasma.pred$v <- influence(age.model)$hat
age.linmodel<-lm(lowplasma ~ age, data = PositivePlasma)
linmodel.pred <- cbind(PositivePlasma,
                   v = influence(age.linmodel)$hat)


(Leverage.plot<-ggplot(PositivePlasma.pred, aes(age, v)) + 
  geom_point() +  
  geom_point(data=linmodel.pred, aes(age, v), color="forestgreen")+
  geom_hline(yintercept = 2*length(age.model$coefficients)/nrow(PositivePlasma), 
             color = "red", size = 1) +
  labs(title = "Leverage vs linear predictor, by Y=0 or Y=1",
       caption = "2(p+1)/n in red, linear & logistic regression leverages in green and black") +
  theme(text = element_text(size = 14)))
#Theyre all 83 years of age.
highv.idx <- which(PositivePlasma.pred$v==max(PositivePlasma.pred$v))
(Q1c.plot<-Q1b.plot+geom_point(data = PositivePlasma.pred[c(highv.idx), ], size = 3, color = "red", shape = 24))
#Question d
PositivePlasma.pred$devres <- influence(age.model)$dev.res
PositivePlasma.pred$devstd <- PositivePlasma.pred$devres/sqrt(1 - PositivePlasma.pred$v)

(devstd.plot<-ggplot(PositivePlasma.pred, aes(age, devstd, color=plasmacat)) + 
    geom_point() +
    geom_point(data = PositivePlasma.pred[c(highv.idx), ], size = 3, color = "red", shape = 24)+
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", size = 1) +
    geom_hline(yintercept = c(-3, 3), linetype = "dotted", size = 1) +
    labs(caption = "\U00B1 2 dashed, \U00B1 3 dotted") +
    theme(text = element_text(size = 14)))

highdevstd.idx <- which(abs(PositivePlasma.pred$devstd)==max(abs(PositivePlasma.pred$devstd)))
(Q1d.plot<-Q1c.plot+geom_point(data = PositivePlasma.pred[c(highdevstd.idx), ], size = 3, color = "blue", shape = 24))
