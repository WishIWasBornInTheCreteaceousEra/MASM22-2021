#Project 3
library(ggplot2)
library(GGally)
library(ggpubr)
library(pROC)
library(multiROC)
library(ResourceSelection)
library(nnet)
library(MASS)

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

#Ordinal
(Qmodel.0<-polr(plasmacat~1, data=PositivePlasma))
(QMax.model <- polr(plasmacat ~ vituse+calories+fiber+alcohol+betadiet+fat+cholesterol+age+sex+smokstat+quetelet, data = PositivePlasma))
#Backwards AIC: 809.1773   ~ vituse + calories + fiber + betadiet + age + sex + smokstat + quetelet
QAICBIntermediate.model<-step(QMax.model, 
                             scope = list(lower = Qmodel.0, upper = QMax.model),
                             direction = "backward")
#Backwards BIC: 814.6367  ~ vituse + betadiet + age + sex + quetelet
QBICBIntermediate.model<-step(QMax.model, 
                             scope = list(lower = Qmodel.0, upper = QMax.model),
                             direction = "backward",
                             k = log(nrow(PositivePlasma)))
#Forward AIC: 808.7619 ~ quetelet + betadiet + vituse + age + smokstat + sex + fiber + fat
QAICFIntermediate.model<-step(Qmodel.0, 
                             scope = list(lower = Qmodel.0, upper = QMax.model),
                             direction = "forward")
#Forwards BIC: 814.6367 ~ quetelet + betadiet + age + vituse + sex
QBICFIntermediate.model<-step(Qmodel.0, 
                             scope = list(lower = Qmodel.0, upper = QMax.model),
                             direction = "forward",
                             k = log(nrow(PositivePlasma)))
# AIC: 809.6859 , 809.488 , 811.6779 
(nosmok.model <- polr(plasmacat ~ quetelet + betadiet + vituse + age + sex + fiber + fat, data = PositivePlasma))
(nofat.model <- polr(plasmacat ~ quetelet + betadiet + vituse + age + smokstat + sex + fiber, data = PositivePlasma))
(nofatnosmok.model <- polr(plasmacat ~ quetelet + betadiet + vituse + age + sex + fiber, data = PositivePlasma))
#Based on our p-values the QAICFIntermediate.model outperforms these three models.
anova(nofatnosmok.model, QAICFIntermediate.model)#p=.0304281
anova(nofatnosmok.model, nosmok.model)#p=.045
anova(nofatnosmok.model, nofat.model)#p=.045
#anova(nofatnosmok.model, Qmodel.0)#p=2.36e-14

#Tests
anova(QBICFIntermediate.model, QAICFIntermediate.model)#AICF Better
anova(QBICBIntermediate.model, QAICFIntermediate.model)#AICF Better
anova(QAICBIntermediate.model, QAICFIntermediate.model)#AICB Better!!


Qmodel.final<-QAICFIntermediate.model
Qsum.final<-summary(Qmodel.final)
Qparam.final<-cbind(beta = Qmodel.final$coefficients, 
      expbeta = exp(Qmodel.final$coefficients),
      exp(confint(Qmodel.final)),
      zeta = Qmodel.final$zeta, 
      expzeta = exp(Qmodel.final$zeta))

y0 <- data.frame(quetelet = rep(seq(16.33114, 50.40333), 2))
y0$vituse <- "Yes, fairly often"
y0$age<- mean(PositivePlasma$age)
y0$fiber<- mean(PositivePlasma$fiber)
y0$betadiet<- mean(PositivePlasma$betadiet)
y0$fat<-mean(PositivePlasma$fat)
y0$smokstat<-"Never"
y0$sex <- c(rep("Male", 50.40333 - 16.33114 + 1), 
            rep("Female", 50.40333 - 16.33114 + 1))

x00<-data.frame(age = rep(seq(19, 83), 2))
x00$vituse <- "Yes, fairly often"
x00$quetelet<- mean(PositivePlasma$quetelet)
x00$fiber<- mean(PositivePlasma$fiber)
x00$betadiet<- mean(PositivePlasma$betadiet)
x00$fat<-mean(PositivePlasma$fat)
x00$smokstat<-"Never"
x00$sex <- c(rep("Male", 83 - 19 + 1), 
            rep("Female", 83 - 19 + 1))

QFinal.pred <- cbind(
  x00,
  predict(Qmodel.final, newdata = x00, type = "prob"),
  yhat = predict(Qmodel.final, newdata = x00))
QQFinal.pred <- cbind(
  y0,
  predict(Qmodel.final, newdata = y0, type = "prob"),
  yhat = predict(Qmodel.final, newdata = y0))

(OrdinalQuetelet.plot<-ggplot(QQFinal.pred, aes(x=quetelet))+
    geom_line(aes(y = High, color = "High"), size = 2) +
    geom_line(aes(y = `Moderately High`, color = "Moderately High"), size = 2) +
    geom_line(aes(y = Low, color = "Low"), size = 2) +
    geom_line(aes(y = `Very Low`, color = "Very Low"), size = 2) +
    expand_limits(y = c(0, 1)) +
    facet_grid(~sex) +
    labs(title = "Ordinal: average patient with varying quetelet and sex",
         y = "probability") +
    theme(text = element_text(size = 14)))

(OrdinalAge.plot<-ggplot(QFinal.pred, aes(x=age))+
    geom_line(aes(y = High, color = "High"), size = 2) +
    geom_line(aes(y = `Moderately High`, color = "Moderately High"), size = 2) +
    geom_line(aes(y = Low, color = "Low"), size = 2) +
    geom_line(aes(y = `Very Low`, color = "Very Low"), size = 2) +
    expand_limits(y = c(0, 1)) +
    facet_grid(~sex) +
    labs(title = "Ordinal: average patient with varying ages and sex",
         y = "probability") +
    theme(text = element_text(size = 14)))


# aic/bic, R2####
# deviance
Qmodel.final$deviance
# total number of parameters (beta and zeta)
Qmodel.final$edf

Qinfo <- cbind(aic = AIC(Qmodel.0, QAICBIntermediate.model, Qmodel.final, QMax.model),
              bic = BIC(Qmodel.0, QAICBIntermediate.model, Qmodel.final, QMax.model),
              R2D = 100*c(0, 
                          1 - QAICBIntermediate.model$deviance/Qmodel.0$deviance,
                          1 - Qmodel.final$deviance/Qmodel.0$deviance,
                          1 - QMax.model$deviance/model.null$deviance),
              R2D.adj = 100*c(0, 
                              1 - (QAICBIntermediate.model$deviance + QAICBIntermediate.model$edf - Qmodel.0$edf)/
                                Qmodel.0$deviance, 
                              1 - (Qmodel.final$deviance + Qmodel.final$edf - Qmodel.0$edf)/
                                Qmodel.0$deviance,
                              1 - (QMax.model$deviance + QMax.model$edf - Qmodel.0$edf)/
                                Qmodel.0$deviance))
round(Qinfo, digits = 1)


#Confusion matrix####

Qpred.final <- cbind(PositivePlasma,
                    yhat = predict(Qmodel.final))
(Qconf.matrix <- table(Qpred.final$plasmacat, Qpred.final$yhat))
table(Qpred.final$plasmacat)
table(Qpred.final$yhat)
sum(Qconf.matrix)

(Qsens <- 100*(diag(Qconf.matrix)/table(Qpred.final$plasmacat)))
(Qprec <- 100*(diag(Qconf.matrix)/table(Qpred.final$yhat)))
(Qacc <- 100*sum(diag(Qconf.matrix)/sum(Qconf.matrix)))



#Multinomial
(model.0 <- multinom(plasmacat ~ 1, data = PositivePlasma))
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
# Sensitivities:
(sens.final <-round(100*diag(conf.final)/table(pred.final$plasmacat), digits = 1))
# precision:
(prec.final <-round(100*diag(conf.final)/table(pred.final$yhat), digits = 1))
# accuracy:
(acc.final<-100*sum(diag(conf.final))/sum(conf.final))


