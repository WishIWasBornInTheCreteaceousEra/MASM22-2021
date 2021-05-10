#Project 2
library(ggplot2)
library(GGally)
library(ggpubr)
library(pROC)
library(ResourceSelection)

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

#Question e:
# Cook's distance####
PositivePlasma.pred$Dcook <- cooks.distance(age.model)

ggplot(PositivePlasma.pred, aes(age, Dcook, color = plasmacat)) +
  geom_point() +
  geom_point(data = PositivePlasma.pred[c(highv.idx), ], size = 3,
            color = "red", shape = 24)+
  geom_point(data = PositivePlasma.pred[c(highdevstd.idx), ], size = 3,
             color = "blue", shape = 24)+
  geom_hline(yintercept = 4/nrow(PositivePlasma), linetype = "dotted",
             size = 1) +
  labs(color = "\U03B2-carotene\nconcentrations", 
       caption = "4/n in black, high leverage highlighted in red and high residual in blue.",
       y="Cook's Distance") +
  theme(text = element_text(size = 14))


##Part 2:##
#Question a:
model.0 <- glm(lowplasma ~ 1, family = "binomial", data = PositivePlasma)
model2.max <- glm(lowplasma ~ age+sex+smokstat+quetelet, family = "binomial", data = PositivePlasma)
background.model<-step(age.model, 
                    scope = list(lower = model.0, upper = model2.max),
                    direction = "both")
(BetaB.lm <- cbind(summary(background.model)$coefficients,ci =confint(background.model)))
(exp(BetaB.lm[,c(1,5,6)]))

PositivePlasma$agecat <- cut(PositivePlasma$age, breaks = c(0, 40, 55, 100))
Background.pred <- cbind(
  PositivePlasma,
  phat = predict(background.model, type = "response"),
  logit = predict(background.model, se.fit = TRUE))

(lowplasmaAGEbg.plot<-ggplot(Background.pred, aes(age, lowplasma)) + 
    geom_point() +
    geom_point(aes(y=phat, color=bmicat))+
    facet_grid(~smokstat)+
    labs(y="probability of having low \U03B2-carotene", 
         caption = "Individual plots seperated by smoking status",
         color = "BMI\nCategories") +
    theme(text = element_text(size = 14)))

(lowplasmaQUETbg.plot<-ggplot(Background.pred, aes(quetelet, lowplasma)) + 
    geom_point() +
    geom_point(aes(y=phat, color=agecat))+
    facet_grid(~smokstat)+
    labs(y="probability of having low \U03B2-carotene", 
         caption = "Individual plots seperated by smoking status",
         color="Age\nGroups") +
    theme(text = element_text(size = 14)))
#Question b:
Background.pred$v <- influence(background.model)$hat
Background.pred$devres <- influence(background.model)$dev.res
Background.pred$devstd <- Background.pred$devres/sqrt(1 - Background.pred$v)
Background.pred$Dcook <- cooks.distance(background.model)


#Leverage Plots
(LeverageBMIbg.plot<-ggplot(Background.pred, aes(age, v, color=bmicat)) + 
    geom_point() +  
    facet_grid(~smokstat)+
    geom_hline(yintercept = 2*length(background.model$coefficients)/nrow(PositivePlasma), 
               color = "red", size = 1) +
    labs(caption = "2(p+1)/n in red",
         color = "BMI\nCategories",
         y="Leverage") +
    theme(text = element_text(size = 14)))

(LeverageAgebg.plot<-ggplot(Background.pred, aes(quetelet, v, color=agecat)) + 
    geom_point() +  
    facet_grid(~smokstat)+
    geom_hline(yintercept = 2*length(background.model$coefficients)/nrow(PositivePlasma), 
               color = "red", size = 1) +
    labs(caption = "2(p+1)/n in red",
         color = "Age\nGroups",
         y="Leverage") +
    theme(text = element_text(size = 14)))

#Deviance plots
(devstdBMIbg.plot<-ggplot(Background.pred, aes(age, devstd, color=bmicat)) + 
    geom_point() +
    geom_point(data = Background.pred[abs(Background.pred$devstd) > 2,], 
                                             color = "red", size = 3)+
    geom_hline(yintercept = 0) + 
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", size = 1) +
    labs(caption = "\U00B1 2 dashed",
         color = "BMI\nCategories",
         y="Standardized Deviance") +
    theme(text = element_text(size = 14)))

(devstdAgebg.plot<-ggplot(Background.pred, aes(quetelet, devstd, color=agecat)) + 
    geom_point() +
    geom_point(data = Background.pred[abs(Background.pred$devstd) > 2,], 
               color = "red", size = 3)+
    geom_hline(yintercept = 0) + 
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", size = 1) +
    labs(caption = "\U00B1 2 dashed",
         color = "Age\nGroups",
         y="Standardized Deviance") +
    theme(text = element_text(size = 14)))

#Cook Plots
(CooksDBMIbg.plot<-ggplot(Background.pred, aes(age, Dcook, color=bmicat)) + 
  geom_point() +
  geom_point(data = Background.pred[abs(Background.pred$devstd) > 2,], 
             color = "blue", size = 3, shape=24)+
  geom_hline(yintercept = 4/nrow(PositivePlasma), linetype = "dotted",
             size = 1) +
  labs(color = "BMI\nCategories", 
       caption = "4/n in black, high residual in blue.",
       y="Cook's Distance") +
  theme(text = element_text(size = 14)))

(CooksDAgebg.plot<-ggplot(Background.pred, aes(quetelet, Dcook, color=agecat)) + 
    geom_point() +
    geom_point(data = Background.pred[abs(Background.pred$devstd) > 2,], 
               color = "blue", size = 3, shape=24)+
    geom_hline(yintercept = 4/nrow(PositivePlasma), linetype = "dotted",
               size = 1) +
    labs(color = "Age\nGroups", 
         caption = "4/n in black, high residual in blue.",
         y="Cook's Distance") +
    theme(text = element_text(size = 14)))
(lowplasmaAGEbg2.plot<-lowplasmaAGEbg.plot+geom_point(data = Background.pred[abs(Background.pred$devstd) > 2,], 
                                                color = "blue", size = 3, shape=24))
(lowplasmaQUETbg2.plot<-lowplasmaQUETbg.plot+geom_point(data = Background.pred[abs(Background.pred$devstd) > 2,], 
                                                   color = "blue", size = 3, shape=24))
#Question c: We use AIC as the number of independent variables is already rather small compared to our degrees of
#freedom thus, limiting it is not necessarily ideal. Furthermore, .backwards is better due to collinearity and the
#fact that the number of samples we have exceed the number of parameters so considering each parameter is beneficial.
model2c.max <- glm(lowplasma ~ vituse+calories+fiber+alcohol+betadiet+fat+cholesterol, family = "binomial", data = PositivePlasma)

#AIC:320.36 with vituse + calories + fiber + betadiet on backwards (See modules file for complete comparison.)
#Explanation above however covers our abses
diet.model<-step(model2c.max, 
                    scope = list(lower = model.0, upper = model2c.max),
                    direction = "backward")

(BetaDiet.lm <- cbind(summary(diet.model)$coefficients,ci =confint(diet.model)))
(exp(BetaDiet.lm[,c(1,5,6)]))

#Testing: We conclude that the diet
#Create object for ease of use.
collect.AIC <- AIC(age.model, background.model, diet.model)
(lnL0 <- logLik(model.0)[1])
(R2CS.max <- 1 - (exp(lnL0))^(2/nrow(PositivePlasma)))
# Collect the log likelihoods L(betahat)
collect.AIC$loglik <- 
  c(logLik(age.model)[1],
    logLik(background.model)[1],
    logLik(diet.model)[1])
# calculate R2_McF;
collect.AIC$R2McF <- 1 - collect.AIC$loglik/lnL0
# calculate R2_McF,adj. Note that p+1 = df (and df.1):
collect.AIC$R2McF.adj <- 1 - (collect.AIC$loglik - (collect.AIC$df - 1)/2)/lnL0
# calculate R2_CS:
collect.AIC$R2CS <- 1 - (exp(lnL0 - collect.AIC$loglik))^(2/nrow(PositivePlasma))
# Calculate R2_N:
collect.AIC$R2N <- collect.AIC$R2CS/R2CS.max

# Show them as % with one decimal value: the p's are 2,5,6. AIC and McFadden's adjusted pseudo R^2 agrees diet is 
#best. Unrelated side note, Cox and Snell's R2 as a theoretical maximum value of less than 1, even for a "perfect"
#model when the outcomes can be categorical. Nagelkerke's R2 fixes that. McFadden's is a comparison to the intercept
#model.
round(100*collect.AIC[, c("R2McF", "R2McF.adj", "R2CS", "R2N")], digits = 1)


#Question d:(Check the modules file for the quick method.)
#We know BIC will punish large models so we use BIC in the backwards configuration to find the lowest state and work
#from there.
model2d.max<- glm(lowplasma ~ vituse+calories+fiber+alcohol+betadiet+fat+cholesterol+age+sex+smokstat+quetelet, family = "binomial", data = PositivePlasma)

#However due to the reasons mentioned before it likely penalizes our model too heavily so we should use AIC in the
#backwards configuration. AIC of 306.9 vituse* + calories + fiber + betadiet + age + smokstat + quetelet
AICintermediate.model<-step(model2d.max, 
                         scope = list(lower = model.0, upper = model2d.max),
                         direction = "backward")
#AIC of 308.91 with vituse* + betadiet + age + quetelet
BICintermediate.model<-step(model2d.max, 
                         scope = list(lower = model.0, upper = model2d.max),
                         direction = "backward",
                         k = log(nrow(PositivePlasma)))
#BIC Test model has an AIC of 306.53 with vituse* + betadiet + age + quetelet + smokstat, *indicates
#vituse of yes often is not used.
BICtest.model<-step(BICintermediate.model, 
                 scope = list(lower = model.0, upper = model2d.max),
                 direction = "both")
#AIC Test model would be identical. Performing a manual analysis we find that fiber, smoking status, and calories 
#are potentially insignificant
Fiberless.model<-glm(lowplasma ~ vituse+calories+betadiet+age+smokstat+quetelet, family = "binomial", data = PositivePlasma)
nosmoke.model<-glm(lowplasma ~ vituse+calories+betadiet+age+fiber+quetelet, family = "binomial", data = PositivePlasma)
nocalories.model<-glm(lowplasma ~ vituse+fiber+betadiet+age+smokstat+quetelet, family = "binomial", data = PositivePlasma)

Nosmokefiber.model<-glm(lowplasma ~ vituse+calories+betadiet+age+quetelet, family = "binomial", data = PositivePlasma)
Nofibercalories.model<-glm(lowplasma ~ vituse+betadiet+age+smokstat+quetelet, family = "binomial", data = PositivePlasma)
Nosmokecalories.model<-glm(lowplasma ~ vituse+fiber+betadiet+age+quetelet, family = "binomial", data = PositivePlasma)

#We care about the adjusted mcfadden pseudo R^2's
AIC.list <- AIC(AICintermediate.model, 
                BICtest.model, 
                Fiberless.model, 
                nosmoke.model,
                nocalories.model, 
                Nosmokefiber.model, 
                Nofibercalories.model, 
                Nosmokecalories.model)
AIC.list$loglikelihoods<-c(logLik(AICintermediate.model)[1],
                  logLik(BICtest.model)[1],
                  logLik(Fiberless.model)[1],
                  logLik(nosmoke.model)[1],
                  logLik(nocalories.model)[1],
                  logLik(Nosmokefiber.model)[1],
                  logLik(Nofibercalories.model)[1],
                  logLik(Nosmokecalories.model)[1])
#We find that our AICintermediate model is the final model, this agrees with the package bestglm.
(AIC.list$R2McF.adj <- 1 - (AIC.list$loglikelihoods/lnL0))
(AIC.list$R2McF.adj <- 1 - (AIC.list$loglikelihoods - (AIC.list$df - 1)/2)/lnL0)

(BetaFinal.lm <- cbind(summary(AICintermediate.model)$coefficients,ci =confint(AICintermediate.model)))
(exp(BetaFinal.lm[,c(1,5,6)]))

#Part 3:
#Question a:
P3.pred<-cbind(PositivePlasma,
               p.age=predict(age.model, type="response"),
               p.background=predict(background.model, type="response"),
               p.diet=predict(diet.model, type="response"),
               p.final=predict(AICintermediate.model, type="response"))


P3.pred$yhat.age <- as.numeric(P3.pred$p.age > 0.5)
P3.pred$yhat.background <- as.numeric(P3.pred$p.background > 0.5)
P3.pred$yhat.diet <- as.numeric(P3.pred$p.diet > 0.5)
P3.pred$yhat.final <- as.numeric(P3.pred$p.final > 0.5)

row.01 <- table(P3.pred$lowplasma)
#Age: Manual entry due to no p<0.5
(col.01.age <- table(P3.pred$yhat.age))
(confusion.age <- table(P3.pred$lowplasma, P3.pred$yhat.age))
(spec.age <- 0 / row.01[1])
(sens.age <- confusion.age[2, 1] / row.01[2])
(accu.age <- confusion.age[2, 1] / sum(confusion.age))
(prec.age <- confusion.age[2, 1] / col.01.age[1])

#Background
(col.01.background <- table(P3.pred$yhat.background))
(confusion.background <- table(P3.pred$lowplasma, P3.pred$yhat.background))
(spec.background <- confusion.background[1, 1] / row.01[1])
(sens.background <- confusion.background[2, 2] / row.01[2])
(accu.background <- sum(diag(confusion.background)) / sum(confusion.background))
(prec.background <- confusion.background[2, 2] / col.01.background[2])

#Diet
(col.01.diet <- table(P3.pred$yhat.diet))
(confusion.diet <- table(P3.pred$lowplasma, P3.pred$yhat.diet))
(spec.diet <- confusion.diet[1, 1] / row.01[1])
(sens.diet <- confusion.diet[2, 2] / row.01[2])
(accu.diet <- sum(diag(confusion.diet)) / sum(confusion.diet))
(prec.diet <- confusion.diet[2, 2] / col.01.diet[2])

#Final
(col.01.final <- table(P3.pred$yhat.final))
(confusion.final <- table(P3.pred$lowplasma, P3.pred$yhat.final))
(spec.final <- confusion.final[1, 1] / row.01[1])
(sens.final <- confusion.final[2, 2] / row.01[2])
(accu.final <- sum(diag(confusion.final)) / sum(confusion.final))
(prec.final <- confusion.final[2, 2] / col.01.final[2])

#Question b
# ROC-curves for all models####
roc.age <- roc(lowplasma ~ p.age, data = P3.pred)
roc.df.age <- coords(roc.age, transpose = FALSE)
roc.df.age$model <- "Age"
roc.background <- roc(lowplasma ~ p.background, data = P3.pred)
roc.df.background <- coords(roc.background, transpose = FALSE)
roc.df.background$model <- "Background"
roc.diet <- roc(lowplasma ~ p.diet, data = P3.pred)
roc.df.diet <- coords(roc.diet, transpose = FALSE)
roc.df.diet$model <- "Diet"
roc.final <- roc(lowplasma ~ p.final, data = P3.pred)
roc.df.final <- coords(roc.final, transpose = FALSE)
roc.df.final$model <- "Final"


roc.df <- rbind(roc.df.age, roc.df.background, roc.df.diet, 
                roc.df.final)

# Plot all the curves, in different colors:
ggplot(roc.df, aes(specificity, sensitivity,
                   color = model)) +
  geom_path(size = 1) +
  coord_fixed() +       # square plotting area
  scale_x_reverse() +   # Reverse scale on the x-axis!
  labs(caption = "ROC-curves for all the models") +
  theme(text = element_text(size = 14))

#Collect AUC and intervals for all the models:
(aucs <- 
    data.frame(
      model = c("Age", "Background", "Diet", "Final"),
      auc = c(auc(roc.age), auc(roc.background), auc(roc.diet), auc(roc.final)),
      lwr = c(ci(roc.age)[1], ci(roc.background)[1],
              ci(roc.diet)[1], ci(roc.final)[1]),
      upr = c(ci(roc.age)[3], ci(roc.background)[3],
              ci(roc.diet)[3], ci(roc.final)[3])))

# Compare the AUC for the models:
roc.test(roc.age, roc.background)
roc.test(roc.age, roc.diet)
roc.test(roc.age, roc.final)
roc.test(roc.background, roc.diet)
roc.test(roc.background, roc.final)
roc.test(roc.diet, roc.final)

#Question c: Look at modules for a simpler more visual method using the cutpointr package, consider using the plots 
#from there to discuss the best threshold points.
MaxSpSe.idx<-c(which.min(sqrt((1-roc.df.age$specificity)^2+(1-roc.df.age$sensitivity)^2)),
               which.min(sqrt((1-roc.df.background$specificity)^2+(1-roc.df.background$sensitivity)^2)),
               which.min(sqrt((1-roc.df.diet$specificity)^2+(1-roc.df.diet$sensitivity)^2)),
               which.min(sqrt((1-roc.df.final$specificity)^2+(1-roc.df.final$sensitivity)^2)))

thresh.age<-roc.df.age$threshold[MaxSpSe.idx[1]]
thresh.background<-roc.df.background$threshold[MaxSpSe.idx[2]]
thresh.diet<-roc.df.diet$threshold[MaxSpSe.idx[3]]
thresh.final<-roc.df.final$threshold[MaxSpSe.idx[4]]

P3.pred$ydhat.age <- as.numeric(P3.pred$p.age > thresh.age)
P3.pred$ydhat.background <- as.numeric(P3.pred$p.background > thresh.background)
P3.pred$ydhat.diet <- as.numeric(P3.pred$p.diet > thresh.diet)
P3.pred$ydhat.final <- as.numeric(P3.pred$p.final > thresh.final)

#Age: Manual entry due to no p<0.5
(col.02.age <- table(P3.pred$ydhat.age))
(confusiond.age <- table(P3.pred$lowplasma, P3.pred$ydhat.age))
(specd.age <- confusiond.age[1, 1] / row.01[1])
(sensd.age <- confusiond.age[2, 2] / row.01[2])
(accud.age <- sum(diag(confusiond.age)) / sum(confusiond.age))
(precd.age <- confusiond.age[2, 2] / col.02.age[2])

#Background
(col.02.background <- table(P3.pred$ydhat.background))
(confusiond.background <- table(P3.pred$lowplasma, P3.pred$ydhat.background))
(specd.background <- confusiond.background[1, 1] / row.01[1])
(sensd.background <- confusiond.background[2, 2] / row.01[2])
(accud.background <- sum(diag(confusiond.background)) / sum(confusiond.background))
(precd.background <- confusiond.background[2, 2] / col.02.background[2])

#Diet
(col.02.diet <- table(P3.pred$ydhat.diet))
(confusiond.diet <- table(P3.pred$lowplasma, P3.pred$ydhat.diet))
(specd.diet <- confusiond.diet[1, 1] / row.01[1])
(sensd.diet <- confusiond.diet[2, 2] / row.01[2])
(accud.diet <- sum(diag(confusiond.diet)) / sum(confusiond.diet))
(precd.diet <- confusiond.diet[2, 2] / col.02.diet[2])

#Final
(col.02.final <- table(P3.pred$ydhat.final))
(confusiond.final <- table(P3.pred$lowplasma, P3.pred$ydhat.final))
(specd.final <- confusiond.final[1, 1] / row.01[1])
(sensd.final <- confusiond.final[2, 2] / row.01[2])
(accud.final <- sum(diag(confusiond.final)) / sum(confusiond.final))
(precd.final <- confusiond.final[2, 2] / col.02.final[2])

#Question d:
# HL using hoslem.test####
Model.list<-list(age.model, background.model, diet.model, AICintermediate.model)
Model.names<-c("Age", "Background", "Dietary", "Final" )
for(i in 1:4){
  exp<-20 #Completely arbitrary number >5
  git<-length(Model.list[[i]]$coefficients)+1
  P3.sort<-P3.pred[order(P3.pred[I(16+i)]), ]
  while(exp>5){
    HL<-hoslem.test(P3.pred$lowplasma, P3.sort[, I(16+i)], g = git)
    exp=min(HL$expected)
    chisq=qchisq(1 - 0.05, git - 2)
    print(sprintf("Smallest expected number in a group of %s model for g=%i is %f", Model.names[i], git, exp))
    if(chisq<HL[1]){
      print(sprintf("The Chi-sq value at a significance of .05 is %f, and chi-sq of the HL test is %f so we REJECT", chisq, HL[[1]]))
    }
    else{
      print(sprintf("The Chi-sq value at a significance of .05 is %f, and chi-sq of the HL test is %f so we ACCEPT", chisq, HL[[1]]))
    }

    HL.df <- data.frame(group = seq(1, git),
                        Obs0 = HL$observed[, 1],
                        Obs1 = HL$observed[, 2],
                        Exp0 = HL$expected[, 1],
                        Exp1 = HL$expected[, 2])
    
    ggplot(HL.df, aes(x = group)) +
      geom_line(aes(y = Obs0, linetype = "observed", color = "Y = 0"), size = 1) +
      geom_line(aes(y = Obs1, linetype = "observed", color = "Y = 1"), size = 1) +
      geom_line(aes(y = Exp0, linetype = "expected", color = "Y = 0"), size = 1) +
      geom_line(aes(y = Exp1, linetype = "expected", color = "Y = 1"), size = 1) +
      labs(title = sprintf("%s Model: Observed and expected in each group", Model.names[i]),
           y = "number of observations") +
      scale_x_continuous(breaks = seq(1, 11)) +
      theme(text = element_text(size = 14))
    #ggsave(sprintf("HLT %s%i.svg", Model.names[i], git))
    git=git+1
  }
}
