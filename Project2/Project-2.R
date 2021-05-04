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

glm(lowplasma ~ smokstat, family = "binomial", data = PositivePlasma)
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
glm(lowplasma ~ smokstat, family = "binomial", data = PositivePlasma)
(model.2 <- glm(lowplasma ~ age, family = "binomial", data = PositivePlasma))

PositivePlasma.pred <- cbind(
  PositivePlasma,
  phat = predict(model.2, type = "response"))

ggplot(PositivePlasma.pred, aes(age, lowplasma)) +
  geom_point() +
  geom_smooth(se = FALSE, linetype = "dashed") +
  geom_line(aes(y = phat), color = "red", size = 1) +
  xlab("age") +
  ylab("low beta-carotene") +
  labs(title = "High PM10 (=1) or Not high PM10 (=0) vs number of cars",
       caption = "red = fitted line, blue dashed = moving average") +
  theme(text = element_text(size = 14))

# beta: log-odds(ratio) with c.i.:
model.2$coefficients
(ci.beta <- confint(model.2))

# Odds (exp(beta0)) and OR, odds ratio, exp(beta1)
exp(model.2$coefficients)
(ci.or <- exp(ci.beta))
#Wald test
summary(model.2)$coefficients
# Since |-0.02835726| < lambda_0.025 = 1.96 we can reject
# H0: beta_1 = 0
#The age has a significant impact on the beta-carotene levels



















