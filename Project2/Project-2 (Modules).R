library(ggplot2)
library(GGally)
library(ggpubr)
library(bestglm)


#Question c: We use AIC as the number of independent variables is already rather small compared to our degrees of
#freedom thus, limiting it is not necessarily ideal. Furthermore, .backwards is better due to collinearity and the
#fact that the number of samples we have exceed the number of parameters so considering each parameter is beneficial.
midDiet.model <- glm(lowplasma ~ vituse+calories+fiber+betadiet, family = "binomial", data = PositivePlasma)
model2c.max <- glm(lowplasma ~ vituse+calories+fiber+alcohol+betadiet+fat+cholesterol, family = "binomial", data = PositivePlasma)

#AIC:320.36 with vituse + calories + fiber + betadiet on backwards (Fastest performance for best result step: 7 or 4 
#computations depending on the starting model. Whereas backwards automatically is 4 computations. )
StepAIC.model<-step(model2c.max, 
                    scope = list(lower = model.0, upper = model2c.max),
                    direction = "backward")

#BIC:342.82 with vituse + betadiet on forward 342.85 on backwards with vituse + calories + fiber + betadiet
#Step yields the same 2 models depending on the initial model and using an intermediate model based off of
#the AIC model yields the a BIC score of 342.85 while reminaing unchanged.
FStepBIC.model<-step(model.0, 
                     scope = list(lower = model.0, upper = model2c.max),
                     direction = "forward",
                     k = log(nrow(PositivePlasma)))
BStepBIC.model<-step(model2c.max, 
                     scope = list(lower = model.0, upper = model2c.max),
                     direction = "backward",
                     k = log(nrow(PositivePlasma)))
#Testing: This should be deleted before we submit, as our explanation for why we choose backwards AIC covers our 
#bases.
bic <- BIC(model.0, age.model, background.model, StepAIC.model, FStepBIC.model, BStepBIC.model)
aic <- AIC(model.0, age.model, background.model, StepAIC.model, FStepBIC.model, BStepBIC.model)
(collect.AIC <- data.frame(aic, bic))

(lnL0 <- logLik(model.0)[1])
(R2CS.max <- 1 - (exp(lnL0))^(2/nrow(PositivePlasma)))
# Collect the log likelihoods L(betahat)
collect.AIC$loglik <- 
  c(logLik(model.0)[1],
    logLik(age.model)[1],
    logLik(background.model)[1],
    logLik(StepAIC.model)[1],
    logLik(FStepBIC.model)[1],
    logLik(BStepBIC.model)[1])
# calculate R2_McF;
collect.AIC$R2McF <- 1 - collect.AIC$loglik/lnL0
# calculate R2_McF,adj. Note that p+1 = df (and df.1):
collect.AIC$R2McF.adj <- 1 - (collect.AIC$loglik - (collect.AIC$df - 1)/2)/lnL0
# calculate R2_CS:
collect.AIC$R2CS <- 1 - (exp(lnL0 - collect.AIC$loglik))^(2/nrow(PositivePlasma))
# Calculate R2_N:
collect.AIC$R2N <- collect.AIC$R2CS/R2CS.max

# Show them as % with one decimal value:
round(100*collect.AIC[, c("R2McF", "R2McF.adj", "R2CS", "R2N")], digits = 1)

#Question d:
#Using a package to see if what we got makes sense
PositivePlasma.Xy = PositivePlasma[,c(1,2,3,4,6,7,8,9,10,11,12,14)]
bglm.AIC = bestglm(Xy = PositivePlasma.Xy, family = binomial, IC = "AIC", 
                   method = "exhaustive")
ans<-bglm.AIC$Subsets
ans$df<-0.5*(ans$AIC+2*ans$logLikelihood)+1
# calculate R2_McF,adj. We find that the largest R2 occurs at row 7 that is 6 variables; note that this isnt the
#maximum AIC. Based on R^2adj,mcf our model is: (lowplasma~age,smokstat,quetelet,vituse,calories,fiber,betadiet)
ans$R2McF.adj <- 1 - (ans$logLikelihood - (ans$df - 1)/2)/lnL0
