model.final<-AICBIntermediate.model
anova(update(model.final, . ~ . - vituse), model.final)
anova(update(model.final, . ~ . - calories), model.final)#3.6%
anova(update(model.final, . ~ . - fiber), model.final)
anova(update(model.final, . ~ . - betadiet), model.final)#2.6%
anova(update(model.final, . ~ . - age), model.final)
anova(update(model.final, . ~ . - sex), model.final)
anova(update(model.final, . ~ . - quetelet), model.final)