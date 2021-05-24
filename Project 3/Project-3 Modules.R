#Ordinal: quetelet + betadiet + vituse + age + smokstat + sex + fiber + fat
Qmodel.final<-QAICFIntermediate.model
anova(update(Qmodel.final, . ~ . - vituse), Qmodel.final)
anova(update(Qmodel.final, . ~ . - smokstat), Qmodel.final)#8.5%
anova(update(Qmodel.final, . ~ . - fiber), Qmodel.final)
anova(update(Qmodel.final, . ~ . - betadiet), Qmodel.final)
anova(update(Qmodel.final, . ~ . - age), Qmodel.final)
anova(update(Qmodel.final, . ~ . - sex), Qmodel.final)
anova(update(Qmodel.final, . ~ . - quetelet), Qmodel.final)
anova(update(Qmodel.final, . ~ . - fat), Qmodel.final)#9.8%

#Mutlinomial
model.final<-AICBIntermediate.model
anova(update(model.final, . ~ . - vituse), model.final)
anova(update(model.final, . ~ . - calories), model.final)#3.6%
anova(update(model.final, . ~ . - fiber), model.final)
anova(update(model.final, . ~ . - betadiet), model.final)#2.6%
anova(update(model.final, . ~ . - age), model.final)
anova(update(model.final, . ~ . - sex), model.final)
anova(update(model.final, . ~ . - quetelet), model.final)


#Multinomial ROC
p.final<-predict(model.final, type = "probs")
pcat<-pred.final$plasmacat
roc.final <- multiclass.roc(pcat ~ p.final)


#Ordinal ROC
Qp.final<-predict(Qmodel.final, type = "probs")
Qpcat<-Qpred.final$plasmacat
Qroc.final <- multiclass.roc(Qpcat ~ Qp.final)

roc.test(Qroc.final, roc.final)
