#Project 3
library(ggplot2)
library(GGally)
library(ggpubr)
library(pROC)
library(ResourceSelection)
library(nnet)

plasma <- read.delim("Data/plasma.txt")
PositivePlasma <- plasma[plasma$betaplasma > 0, ]
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

#Making our categories
(Quants=quantile(PositivePlasma$betaplasma))
PositivePlasma$plasmacat <- cut(PositivePlasma$betaplasma, 
                                breaks = c(Quants[1]-1, Quants[2], Quants[3], Quants[4], Quants[5]+1),
                                labels=c("Very Low", "Low", "Moderately High", "High"))
table(PositivePlasma$plasmacat)

(model.0 <- multinom(plasmacat ~ 1, data = PositivePlasma))
(background.model <- multinom(plasmacat ~ age+sex+smokstat+quetelet, data = PositivePlasma))
(dietary.model <- multinom(plasmacat ~ vituse+calories+fiber+alcohol+betadiet+fat+cholesterol, data = PositivePlasma))
(Max.model <- multinom(plasmacat ~ vituse+calories+fiber+alcohol+betadiet+fat+cholesterol+age+sex+smokstat+quetelet, data = PositivePlasma))

#Backwards AIC: 818.1562  ~ vituse + calories + fiber + betadiet + age + sex + quetelet
AICBIntermediate.model<-step(Max.model, 
                            scope = list(lower = model.0, upper = Max.model),
                            direction = "backward")
#Backwards BIC: 818.1562  ~ vituse + calories + fiber + betadiet + age + sex + quetelet
BICBIntermediate.model<-step(Max.model, 
                            scope = list(lower = model.0, upper = Max.model),
                            direction = "backward",
                            k = log(nrow(PositivePlasma)))
#Forward AIC: 818 ~ vituse + quetelet + age + betadiet + sex + fiber + calories
AICFIntermediate.model<-step(model.0, 
                            scope = list(lower = model.0, upper = Max.model),
                            direction = "forward")
#Forwards BIC: 857 ~vituse
BICFIntermediate.model<-step(model.0, 
                             scope = list(lower = model.0, upper = Max.model),
                             direction = "forward",
                             k = log(nrow(PositivePlasma)))

sum.final<-summary(AICBIntermediate.model)
(beta <- sum.final$coefficients)
(se.beta <- sum.final$standard.errors)
(z.value <- beta/se.beta)
(P.value <- pnorm(abs(z.value), lower.tail = FALSE))

# AIC: 820.6893 , 821.3737 , 823.7628 
(nocalories.model <- multinom(plasmacat ~ vituse + fiber + betadiet + age + sex + quetelet, data = PositivePlasma))
(nobetadiet.model <- multinom(plasmacat ~ vituse + fiber + calories + age + sex + quetelet, data = PositivePlasma))
(nocalnoBD.model <- multinom(plasmacat ~ vituse + fiber + age + sex + quetelet, data = PositivePlasma))


# LR (deviance) tests#### Better than all but Max.model, however we cannot conclude that
#it is worse.
# against no calories model
anova(AICBIntermediate.model, nocalories.model)
# against no beta diet model
anova(AICBIntermediate.model, nobetadiet.model)
# against no calories no betadeit model
anova(AICBIntermediate.model, nocalnoBD.model)
# against intercept model
anova(AICBIntermediate.model, model.0)
#against max model
anova(AICBIntermediate.model, Max.model)

#Compute R2 Scores: Conclude that AIC Intermediate.model is the best.
info <- cbind(aic = AIC(AICBIntermediate.model, nocalories.model, nobetadiet.model, nocalnoBD.model, Max.model, model.0),
              bic = BIC(AICBIntermediate.model, nocalories.model, nobetadiet.model, nocalnoBD.model, Max.model, model.0))
info$r2 <- round(100*c(
  1 - AICBIntermediate.model$deviance/model.0$deviance,
  1 - nocalories.model$deviance/model.0$deviance,
  1 - nobetadiet.model$deviance/model.0$deviance,
  1 - nocalnoBD.model$deviance/model.0$deviance,
  1 - Max.model$deviance/model.0$deviance,
  0), digits = 1)
info$r2.adj <- round(100*c(
  1 - (AICBIntermediate.model$deviance + (AICBIntermediate.model$edf - model.0$edf))/model.0$deviance,
  1 - (nocalories.model$deviance + (nocalories.model$edf - model.0$edf))/model.0$deviance,
  1 - (nobetadiet.model$deviance + (nobetadiet.model$edf - model.0$edf))/model.0$deviance,
  1 - (nocalnoBD.model$deviance + (nocalnoBD.model$edf - model.0$edf))/model.0$deviance,
  1 - (Max.model$deviance + (Max.model$edf - model.0$edf))/model.0$deviance,
  0),
  digits = 1)


#From here on we are working with our "Final Model"
#odds ratios####
OR <- exp(beta)
# Confidence intervals for OR:
ci <- exp(confint(AICBIntermediate.model))
#Nice Tables for the categories
(Low.dat<-cbind(beta = round(beta["Low", ], digits = 2),
                P.value = round(P.value["Low", ], digits = 3),
                OR = round(OR["Low", ], digits = 2), 
                round(ci[, , "Low"], digits = 2)))
(ModHigh.dat<-cbind(beta = round(beta["Moderately High", ], digits = 2),
               P.value = round(P.value["Moderately High", ], digits = 3),
               OR = round(OR["Moderately High", ], digits = 2), 
               round(ci[, , "Moderately High"], digits = 2)))
(High.dat<-cbind(beta = round(beta["High", ], digits = 2),
               P.value = round(P.value["High", ], digits = 3),
               OR = round(OR["High", ], digits = 2), 
               round(ci[, , "High"], digits = 2)))

x0 <- data.frame(age = rep(seq(19, 83), 2))
# the same ses, socst and science for all:
x0$vituse <- "Yes, fairly often"
x0$calories<- mean(PositivePlasma$calories)
x0$quetelet<- mean(PositivePlasma$quetelet)
x0$fiber<- mean(PositivePlasma$fiber)
x0$betadiet<- mean(PositivePlasma$betadiet)

x0$sex <- c(rep("Male", 83 - 19 + 1), 
               rep("Female", 83 - 19 + 1))
Final.pred <- cbind(x0,
                 predict(AICBIntermediate.model, newdata=x0, type = "probs"))

(something.plot<-ggplot(Final.pred, aes(x=age))+
  geom_line(aes(y = High, color = "High"), size = 2) +
  geom_line(aes(y = `Moderately High`, color = "Moderately High"), size = 2) +
  geom_line(aes(y = Low, color = "Low"), size = 2) +
  geom_line(aes(y = `Very Low`, color = "Very Low"), size = 2) +
  expand_limits(y = c(0, 1)) +
  facet_grid(~sex) +
  labs(title = "Multinomial: average patient with varying ages and sex",
       y = "probability") +
  theme(text = element_text(size = 14)))

(model.final <- multinom(plasmacat ~ vituse+calories+fiber+betadiet+age+sex+quetelet, data = PositivePlasma))
pred.final <- cbind(PositivePlasma,
                    predict(model.final, type = "probs"),
                    yhat = predict(model.final))
# confusion matrix####
table(pred.final$plasmacat)
table(pred.final$yhat)

(conf.final <- table(pred.final$plasmacat, pred.final$yhat))
# Sensitivities:THIS NEEDS TO BE FIXED.
(sens.final <-round(100*diag(conf.final)/table(pred.final$plasmacat), digits = 1))
# Specificity:
(spec.final <-round(100*diag(conf.final)/table(pred.final$plasmacat), digits = 1))
# precision:
(prec.final <-round(100*diag(conf.final)/table(pred.final$yhat), digits = 1))
# accuracy:
(acc.final<-100*sum(diag(conf.final))/sum(conf.final))


