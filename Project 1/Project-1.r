#Project 1
library(ggplot2)

plasma <- read.delim("Data/plasma.txt")
head(plasma)
summary(plasma)
# At least one data point was 0, this needs to be taken 
PositivePlasma <- plasma[plasma$betaplasma > 0, ]
summary(PositivePlasma)

#Linear fit of betaplasma as function of age
plasma.LinearModel <- lm(betaplasma ~ age, data = PositivePlasma)
summary(plasma.LinearModel)
confint(plasma.LinearModel, level=0.95)

#Linear fit of log(betaplasma) as function of age
plasma.LogModel <- lm(log(betaplasma) ~ age, data = PositivePlasma)
summary(plasma.LogModel)
confint(plasma.LogModel, level=0.95)

#Create fit object for the two model
y0.lin<-cbind(PositivePlasma,
              fit=predict(plasma.LinearModel),
              conf=predict(plasma.LinearModel, interval="confidence"),
              pred=predict(plasma.LinearModel, interval="prediction"))

y0.log<-cbind(PositivePlasma,
              fit=predict(plasma.LogModel),
              conf=predict(plasma.LogModel, interval="confidence"),
              pred=predict(plasma.LogModel, interval="prediction"))

#Plots for both betaplasma varying as a function of age, and log(betaplasma) varying as a function of age:
(ggplot(data = y0.lin, aes(x = age, y = betaplasma)) +
  geom_point(size = 3) +   xlab("age [year]") +
  ylab("plasma beta-carotene [ng/ml]")+
  theme(text = element_text(size = 20))+
  geom_line(aes(y=fit))+
    geom_ribbon(aes(ymin = conf.lwr, ymax = conf.upr), alpha = 0.2)+
  geom_line(aes(y = pred.lwr), color = "red", linetype = 2)+
  geom_line(aes(y = pred.upr), color = "red", linetype = 2)
)

(ggplot(data = y0.log, aes(x = age, y = log(betaplasma)))+
  geom_point(size = 3)+
  xlab("age [year]") +
  ylab("log(plasma beta-carotene) [log(ng/ml)]")+
  theme(text = element_text(size = 20))+
  geom_line(aes(y=fit))+
  geom_ribbon(aes(ymin = conf.lwr, ymax = conf.upr), alpha = 0.2)+
  geom_line(aes(y = pred.lwr), color = "red", linetype = 2)+
  geom_line(aes(y = pred.upr), color = "red", linetype = 2)
)


#QQ & residual plots linear model
y0.lin$e <- plasma.LinearModel$residuals

ggplot(y0.lin, aes(x = fit, y = e))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "loess")

ggplot(y0.lin,aes(sample = e))+
  geom_qq()+ 
  geom_qq_line()

ggplot(y0.lin,aes(x = e))+ 
  geom_histogram(bins = 20)

#QQ & residual plots Log model
y0.log$e <- plasma.LogModel$residuals

ggplot(y0.log, aes(x = fit, y = e))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "loess")

ggplot(y0.log,aes(sample = e))+
  geom_qq()+ 
  geom_qq_line()

ggplot(y0.log,aes(x = e))+ 
  geom_histogram(bins = 20)


#Prediction intervals based on age
x0<-data.frame(age = 30)
confint0<-cbind(x0,round(predict(plasma.LogModel, newdata = x0, interval = "confidence"),digits = 12))
predint0<-cbind(x0,round(predict(plasma.LogModel, newdata = x0, interval = "prediction"),digits = 12))

x1<-data.frame(age = 70)
confint1<-cbind(x1,round(predict(plasma.LogModel, newdata = x1, interval = "confidence"),digits = 12))
predint1<-cbind(x1,round(predict(plasma.LogModel, newdata = x1, interval = "prediction"),digits = 12))

x2<-data.frame(age = 31)
confint2<-cbind(x2,round(predict(plasma.LogModel, newdata = x2, interval = "confidence"),digits = 12))
predint2<-cbind(x2,round(predict(plasma.LogModel, newdata = x2, interval = "prediction"),digits = 12))

x3<-data.frame(age = 71)
confint3<-cbind(x3,round(predict(plasma.LogModel, newdata = x3, interval = "confidence"),digits = 12))
predint3<-cbind(x3,round(predict(plasma.LogModel, newdata = x3, interval = "prediction"),digits = 12))

#Part 2, Frequency tables, need to turn categorigal variables into factors first
(
  PositivePlasma$sex <- factor(PositivePlasma$sex, 
                               levels = c(1, 2),
                               labels = c("Male", "Female"))
)
(
  PositivePlasma$smokstat <- factor(PositivePlasma$smokstat,
                           levels = c(1, 2, 3),
                           labels = c("Never", "Former", "Current"))
)
(
  PositivePlasma$bmicat <- factor(PositivePlasma$bmicat,
                         levels = c(1,2,3,4),
                         labels = c("Underweight", "Normal", "Overweight", "Obese"))
)

(frequency <- table(PositivePlasma$smokstat,PositivePlasma$bmicat,PositivePlasma$sex))
(smokefreq <- margin.table(frequency,1))
(bmifreq <- margin.table(frequency,2))
(sexfreq <- margin.table(frequency,3))

#2b) model where beta-carotine depends only on bmi. Show how releveling to Normal matters.
(
bmi.model<- lm(betaplasma~bmicat, data=PositivePlasma)
)
summary(bmi.model)
(x0 <- data.frame(
  x = c(1,2,3,4),
  bmicat = rep(c("Underweight", "Normal", "Overweight", "Obese"),1)))
(y0 <- cbind(
  x0, 
  predict(bmi.model, x0, se.fit = TRUE),
  conf = predict(bmi.model, x0, interval = "confidence"),
  pred = predict(bmi.model, x0, interval = "prediction"))
)
y0$conf.fit <- y0$pred.fit <- NULL
y0$se.pred <- sqrt(y0$se.fit^2 + y0$residual.scale^2)
y0

(
  bmicat.lines <- ggplot(data=PositivePlasma, aes(x=bmicat, y=betaplasma))+
    geom_point(aes(color=bmicat, shape=bmicat), size=3)+
    xlab("BMI Category") +
    ylab("Plasma \U03B2-carotene Levels [ng/ml]")+
    geom_segment(aes(x=x-.15, y=fit, xend=x+.15, yend=fit, color=bmicat), y0, size=1.25)+
    geom_segment(aes(x=x-.15, y=conf.lwr, xend=x+.15, yend=conf.lwr, color=bmicat), y0, linetype="dotdash", size=1.25)+
    geom_segment(aes(x=x-.15, y=conf.upr, xend=x+.15, yend=conf.upr, color=bmicat), y0, linetype="dotdash", size=1.25)+
    geom_segment(aes(x=x-.15, y=pred.lwr, xend=x+.15, yend=pred.lwr, color=bmicat), y0, linetype="dashed", size=1.25)+
    geom_segment(aes(x=x-.15, y=pred.upr, xend=x+.15, yend=pred.upr, color=bmicat), y0, linetype="dashed", size=1.25)
)

#releveling to "Normal" because it has most candidates
PositivePlasma$bmicat <- relevel(PositivePlasma$bmicat,"Normal")
bmimodel_relevel <- lm(betaplasma~bmicat,data=PositivePlasma)
(
  bmi.model_relevel <- lm(betaplasma~bmicat, data=PositivePlasma)
)
summary(bmi.model_relevel)
(x1 <- data.frame(
  x = c(1,2,3,4),
  bmicat = rep(c("Underweight", "Normal", "Overweight", "Obese"),1)))
(y1 <- cbind(
  x1, 
  predict(bmi.model_relevel, x1, se.fit = TRUE),
  conf = predict(bmi.model_relevel, x1, interval = "confidence"),
  pred = predict(bmi.model_relevel, x1, interval = "prediction"))
)
y1$conf.fit <- y1$pred.fit <- NULL
y1$se.pred <- sqrt(y1$se.fit^2 + y1$residual.scale^2)
y1

(
  bmicat.lines <- ggplot(data=PositivePlasma, aes(x=bmicat, y=betaplasma))+
    geom_point(aes(color=bmicat, shape=bmicat), size=3)+
    xlab("BMI Category") +
    ylab("Plasma \U03B2-carotene Levels [ng/ml]")+
    geom_segment(aes(x=x-.15, y=fit, xend=x+.15, yend=fit, color=bmicat), y1, size=1.25)+
    geom_segment(aes(x=x-.15, y=conf.lwr, xend=x+.15, yend=conf.lwr, color=bmicat), y1, linetype="dotdash", size=1.25)+
    geom_segment(aes(x=x-.15, y=conf.upr, xend=x+.15, yend=conf.upr, color=bmicat), y1, linetype="dotdash", size=1.25)+
    geom_segment(aes(x=x-.15, y=pred.lwr, xend=x+.15, yend=pred.lwr, color=bmicat), y1, linetype="dashed", size=1.25)+
    geom_segment(aes(x=x-.15, y=pred.upr, xend=x+.15, yend=pred.upr, color=bmicat), y1, linetype="dashed", size=1.25)
)


#Question 2C
PositivePlasma$sex<-relevel(PositivePlasma$sex, "Female")
PositivePlasma$smokstat<-relevel(PositivePlasma$smokstat, "Never")
Q2c.model<- lm(betaplasma~age+bmicat+smokstat+sex, data=PositivePlasma)
BetaRG.lm <- cbind(summary(Q2c.model)$coefficients,ci =confint(Q2c.model))

#Q2c.1 & 4 get p-value=3.399e-05<.05 and thus, we reject the null hypothesis that the additional variables could be 0. For 4, We accept the null hypothesis that underweight is insignificant given all variables. 
summary(Q2c.model)
#Q2c.2 ANOVA for age, Our new model is better as the p-value is less than .05.
plasma.anova <- anova(plasma.LinearModel, Q2c.model)
(Fvalue<-plasma.anova$F[2])
(Pvalue_RG<-plasma.anova$"Pr(>F)"[2])
#Q2c.3 We can reject or accept the value the variables bring by their p-value. bmicat and smokstat matter but sex and age dont.
anova(Q2c.model)

#sup
