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
bmi.model<- lm(log(betaplasma)~bmicat, data=PositivePlasma)
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
  bmicat.lines <- ggplot(data=PositivePlasma, aes(x=bmicat, y=log(betaplasma)))+
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
(
  bmi.model_relevel <- lm(log(betaplasma)~bmicat, data=PositivePlasma)
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
  bmicat_relevel.lines <- ggplot(data=PositivePlasma, aes(x=bmicat, y=log(betaplasma)))+
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
Q2c.model<- lm(log(betaplasma)~age+bmicat+smokstat+sex, data=PositivePlasma)
BetaRG.lm <- cbind(summary(Q2c.model)$coefficients,ci =confint(Q2c.model))

#Q2c.1 & 4 get p-value=3.399e-05<.05 and thus, we reject the null hypothesis that the additional variables could be 0. For 4, We accept the null hypothesis that underweight is insignificant given all variables. 
summary(Q2c.model)
#Q2c.2 ANOVA for age, Our new model is better as the p-value is less than .05.
plasma.anova <- anova(plasma.LogModel, Q2c.model)
(Fvalue<-plasma.anova$F[2])
(Pvalue_RG<-plasma.anova$"Pr(>F)"[2])
#Q2c.3 We can reject or accept the value the variables bring by their p-value. bmicat and smokstat matter but sex and age dont.
anova(Q2c.model)

#Q2d:
Full.model<-cbind(PositivePlasma,
                        fit=predict(Q2c.model),
                        conf=predict(Q2c.model, interval="confidence"),
                        pred=predict(Q2c.model, interval="prediction"))
(
  ggplot(data = Full.model, aes(x = age, y = log(betaplasma), color = sex))+
    geom_point(size = 3)+
    xlab("age [year]") +
    ylab("log(plasma beta-carotene) [log(ng/ml)]")+
    theme(text = element_text(size = 20))+
    geom_line(aes(y=fit))+
    geom_ribbon(aes(ymin = conf.lwr, ymax = conf.upr), alpha = 0.2)+
    geom_line(aes(y = pred.lwr), linetype = 2)+
    geom_line(aes(y = pred.upr), linetype = 2)+
    facet_grid(smokstat ~ relevel(bmicat, "Underweight"))
)

# The average age for males is 62:
MaleAverageAge <- round(mean(PositivePlasma$age[PositivePlasma$sex=="Male"]))
#Prediction intervals based on age: not sure yet whether can be trusted. "Small" intervals but is this enough?
x2d<-data.frame(sex = "Male", age = MaleAverageAge, smokstat = "Former", bmicat = "Underweight")
(confint2d<-cbind(x2d,round(predict(Q2c.model, newdata = x2d, interval = "confidence"),digits = 12)))
(predint2d<-cbind(x2d,round(predict(Q2c.model, newdata = x2d, interval = "prediction"),digits = 12)))


#Q2e
Q2e.model<- lm(log(betaplasma)~age+quetelet+smokstat+sex, data=PositivePlasma)
BetaQuet.coef <- cbind(summary(Q2e.model)$coefficients,ci =confint(Q2e.model))
##Beta's and exponentiated betas: There is no significant change in the other variables' parameters
(BetaQuet.coef[, c(1, 5,6)])
(exp((BetaQuet.coef[, c(1, 5,6)])))
##Predictions for man and woman 30 yr old, former smokers, normal bmi of 20
F30B20Disc<-data.frame(sex = "Female", age = 30, smokstat = "Former", bmicat = "Normal")
(confint2eDF<-cbind(F30B20Disc,round(predict(Q2c.model, newdata = F30B20Disc, interval = "confidence"),digits = 12)))

M30B20Disc<-data.frame(sex = "Male", age = 30, smokstat = "Former", bmicat = "Normal")
(confint2eDM<-cbind(M30B20Disc,round(predict(Q2c.model, newdata = M30B20Disc, interval = "confidence"),digits = 12)))

F30B20Cont<-data.frame(sex = "Female", age = 30, smokstat = "Former", quetelet = 20)
(confint2eCF<-cbind(F30B20Cont,round(predict(Q2e.model, newdata = F30B20Cont, interval = "confidence"),digits = 12)))

M30B20Cont<-data.frame(sex = "Male", age = 30, smokstat = "Former", quetelet = 20)
(confint2eCM<-cbind(M30B20Cont,round(predict(Q2e.model, newdata = M30B20Cont, interval = "confidence"),digits = 12)))
##Obese (32) Same people
F30B32Disc<-data.frame(sex = "Female", age = 30, smokstat = "Former", bmicat = "Obese")
(confint2eDFO<-cbind(F30B32Disc,round(predict(Q2c.model, newdata = F30B32Disc, interval = "confidence"),digits = 12)))

M30B32Disc<-data.frame(sex = "Male", age = 30, smokstat = "Former", bmicat = "Obese")
(confint2eDMO<-cbind(M30B32Disc,round(predict(Q2c.model, newdata = M30B32Disc, interval = "confidence"),digits = 12)))

F30B32Cont<-data.frame(sex = "Female", age = 30, smokstat = "Former", quetelet = 32)
(confint2eCFO<-cbind(F30B32Cont,round(predict(Q2e.model, newdata = F30B32Cont, interval = "confidence"),digits = 12)))

M30B32Cont<-data.frame(sex = "Male", age = 30, smokstat = "Former", quetelet = 32)
(confint2eCMO<-cbind(M30B32Cont,round(predict(Q2e.model, newdata = M30B32Cont, interval = "confidence"),digits = 12)))
#Relative difference of Obese and normal+confint
(exp(abs(confint2eDFO[5:7]-confint2eDF[5:7]))/exp(confint2eDF[5:7]))#Female BMI
(exp(abs(confint2eDM[5:7]-confint2eDMO[5:7]))/exp(confint2eDM[5:7]))#Male BMI
(exp(abs(confint2eCF[5:7]-confint2eCFO[5:7]))/exp(confint2eCF[5:7]))#Female Quetelet
(exp(abs(confint2eCM[5:7]-confint2eCMO[5:7]))/exp(confint2eCM[5:7]))#Male Quetelet
#Q2f
Q2f.model<- lm(log(betaplasma)~age+quetelet+bmicat+smokstat+sex, data=PositivePlasma)
(summary(Q2f.model))


#Part 3
#Question a:
sum.0<-summary(Q2c.model)
sum.1<-summary(Q2e.model)

#Larger is better
(collect.R2s <- data.frame(
  nr = seq(1, 2),
  model = c("BMIcat", "Quetelet"),
  R2 = c(sum.0$r.squared,
         sum.1$r.squared),
  R2.adj = c(sum.0$adj.r.squared,
             sum.1$adj.r.squared)))

#Smaller is better
(collect.AIC <- data.frame(
  nr = seq(1, 2),
  model = c("BMIcat", "Quetelet"),
  AIC(Q2c.model, Q2e.model),
  BIC(Q2c.model, Q2e.model)))

#Question b:
contx <- plasma[, c("age", "quetelet", "calories", "fat", "fiber",
                    "alcohol", "cholesterol", "betadiet")]
#(Fat, calories) & (Fat, Cholesterol) highly correlated. Potential issue with alcohol and other variables (flat)
ggpairs(contx)

(
  ggplot(data = PositivePlasma, aes(x = fat, y = calories))+
    geom_point(size = 3)+
    xlab("Fat Intake [g/day]") +
    ylab("Cholesterol Intake [mg/day]")+
    theme(text = element_text(size = 20))
)

(
  ggplot(data = PositivePlasma, aes(x = fat, y = cholesterol))+
    geom_point(size = 3)+
    xlab("Fat Intake [g/day]") +
    ylab("Caloric Intake [Cal/day]")+
    theme(text = element_text(size = 20))
)
(
PositivePlasma$vituse <- factor(PositivePlasma$vituse,
                                    levels = c(1, 2, 3),
                                    labels = c("Yes, fairly often", "Yes, not often", "No")))
#Use sex as well to show break down into male and female so that we can choose between Yes, often and No, since they have similar numbers
#We want to choose female max (I think)
(vitusefreq <- table(PositivePlasma$vituse, PositivePlasma$sex))
PositivePlasma$vituse <- relevel(PositivePlasma$vituse,"Yes, fairly often")

#Question c:
Q3c.model<- lm(betaplasma~vituse+calories+fat+fiber+alcohol+cholesterol+betadiet, data=PositivePlasma)
Q3c.pred <- cbind(
  PositivePlasma, 
  fit = predict(Q3c.model),
  r = rstudent(Q3c.model))
Q3c.pred$v <- influence(Q3c.model)$hat

HiV<-2*length(Q3c.model$coefficients)/nrow(PositivePlasma)
Alcoholic<-which(PositivePlasma$alcohol==max(PositivePlasma$alcohol))
#Number of non-drinkers
(length(which(PositivePlasma$alcohol==0)))

p1<-ggplot(Q3c.pred, aes(x = vituse, y = v)) +
  geom_point(size = 2)  +
  geom_point(data = Q3c.pred[abs(Q3c.pred$v) > HiV,], 
             color = "red", size = 3)+
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  geom_hline(yintercept = 1/nrow(PositivePlasma)) +
  geom_hline(yintercept = HiV, 
             color = "red") +
  theme(axis.title.y=element_blank())+
  labs(x = "Vitamin Usage") +
  theme(text = element_text(size = 18))

p2<-ggplot(Q3c.pred, aes(x = calories, y = v)) +
  geom_point(size = 2)  +
  geom_point(data = Q3c.pred[abs(Q3c.pred$v) > HiV,], 
             color = "red", size = 3)+
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  geom_hline(yintercept = 1/nrow(PositivePlasma)) +
  geom_hline(yintercept = HiV, 
             color = "red") +
  theme(axis.title.y=element_blank())+
  labs(x = "Caloric Intake [Cal/day]") +
  theme(text = element_text(size = 18))

p3<-ggplot(Q3c.pred, aes(x = fat, y = v)) +
  geom_point(size = 2)  +
  geom_point(data = Q3c.pred[abs(Q3c.pred$v) > HiV,], 
             color = "red", size = 3)+
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  geom_hline(yintercept = 1/nrow(PositivePlasma)) +
  geom_hline(yintercept = HiV, 
             color = "red") +
  theme(axis.title.y=element_blank())+
  labs(x = "Fat Intake [mg/day]") +
  theme(text = element_text(size = 18))

p4<-ggplot(Q3c.pred, aes(x = fiber, y = v)) +
  geom_point(size = 2)  +
  geom_point(data = Q3c.pred[abs(Q3c.pred$v) > HiV,], 
             color = "red", size = 3)+
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  geom_hline(yintercept = 1/nrow(PositivePlasma)) +
  geom_hline(yintercept = HiV, 
             color = "red") +
  theme(axis.title.y=element_blank())+
  labs(x = "Dietary Fiber Intake [g/day]") +
  theme(text = element_text(size = 18))

p5<-ggplot(Q3c.pred, aes(x = cholesterol, y = v)) +
  geom_point(size = 2)  +
  geom_point(data = Q3c.pred[abs(Q3c.pred$v) > HiV,], 
             color = "red", size = 3)+
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  geom_hline(yintercept = 1/nrow(PositivePlasma)) +
  geom_hline(yintercept = HiV, 
             color = "red") +
  theme(axis.title.y=element_blank())+
  labs(x = "Cholesterol Intake [mg/day]") +
  theme(text = element_text(size = 18))

p6<-ggplot(Q3c.pred, aes(x = alcohol, y = v)) +
  geom_point(size = 2)  +
  geom_point(data = Q3c.pred[abs(Q3c.pred$v) > HiV,], 
             color = "red", size = 3)+
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  geom_hline(yintercept = 1/nrow(PositivePlasma)) +
  geom_hline(yintercept = HiV, 
             color = "red") +
  theme(axis.title.y=element_blank())+
  labs(x = "Alcohol Intake [drinks/week]") +
  theme(text = element_text(size = 18))

p7<-ggplot(Q3c.pred, aes(x = betadiet, y = v)) +
  geom_point(size = 2)  +
  geom_point(data = Q3c.pred[abs(Q3c.pred$v) > HiV,], 
             color = "red", size = 3)+
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  geom_hline(yintercept = 1/nrow(PositivePlasma)) +
  geom_hline(yintercept = HiV, 
             color = "red") +
  theme(axis.title.y=element_blank())+
  labs(x = "Dietary beta-carotene Intake [\U03BCg/day]") +
  theme(text = element_text(size = 18))


Plots<-ggarrange(p1, p2, p3, p4, p5, p6, p7, labels="AUTO")
annotate_figure(Plots,
                top = text_grob("Leverage of the dietary variables", color = "black", face = "bold", size = 14),
                bottom = text_grob("y = 1/n (black) and 2(p+1)/n (red)", color = "black",
                                   hjust = 1, x = 1, face = "italic", size = 10),
                left = text_grob("Leverage", color = "black", rot = 90),
)
#Taking the log of alcohol consumption is a bad idea as we lose 110 non drinkers data. That said, we generally take 
#the log (or any transform) to normalize the distribution of the variable and stabilize its variance. We can see
#that the alcohol intake is indeed non-normally distributed, but its p-value is .222 and thus, rather irrelevant.
#Consequently, eliminating this individual from the dataset is likely a better recourse. (Rawlings)

#Question d: Clearly,from the plot we can see that the variance increases as beta-carotene levels increase. Why?
(
  ggplot(Q3c.pred, aes(x = fit, y = r)) +
    geom_point(size = 3) +
    geom_smooth(color = "red")+
    geom_point(data = Q3c.pred[Alcoholic, ], 
               color = "blue", shape = 24, size = 3)+
    geom_hline(yintercept = c(-2, 0, 2)) +
    geom_hline(yintercept = c(-3, 3), linetype = 2) +
    geom_point(data = Q3c.pred[abs(Q3c.pred$r) > 3,], color = "red", size = 4) +
    xlab("Predicted \U03B2-carotene plasma levels [ng/ml]") +
    ylab("Studentized residual") +
    labs(title = "Studentized residuals vs predicted \U03B2-carotene plasma levels") +
    theme(text = element_text(size = 18))
)
#These two aren't necessary but are nice to see distribution and fit a bit better. 
(
  ggplot(Q3c.pred,aes(sample = r))+
    geom_qq()+ 
    geom_qq_line()+
    labs(title="Q-Q plot of dietary model")+
    theme(text = element_text(size = 18))
)
(
  ggplot(Q3c.pred,aes(x = r))+
    labs(title="Distribution of residuals in the dietary model", x="Studentized residual")+
    theme(text = element_text(size = 18))+
    geom_histogram(bins = 30)
)

ExResid<-which(abs(Q3c.pred$r) > 3)
PositivePlasma[c(ExResid),]

