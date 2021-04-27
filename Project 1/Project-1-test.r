library(tidyverse)
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)



iris.gathered<-pivot_longer(PositivePlasma[,-c(1,2,3,4,5)], cols = -vituse, names_to = "var", values_to = "val")
head(iris.gathered, 3)
(
  ggplot(iris.gathered, aes(x = val, y = vituse)) +
    geom_point(size = 2)  +
    geom_point(data = Q3c.pred[abs(Q3c.pred$v) > HiV,], 
               color = "red", size = 3)+
    geom_point(data = Q3c.pred[Alcoholic, ], 
               color = "blue", shape = 24, size = 3) +
    geom_hline(yintercept = 1/nrow(PositivePlasma)) +
    geom_hline(yintercept = HiV, 
               color = "red") +
    facet_wrap(~var, scale="free")+
    theme(axis.title.y=element_blank())+
    labs(x = "Dietary beta-carotene Intake [\U03BCg/day]") +
    theme(text = element_text(size = 18))+
    scale_color_viridis_d()
)


ggplot(Q3c.pred, aes(x = betadiet, y = v)) +
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