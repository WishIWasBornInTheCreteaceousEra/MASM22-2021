pa1<-ggplot(Q3c.pred, aes(x = vituse, y = alcohol)) +
  geom_point(size = 2)+
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  theme(axis.title.y=element_blank())+
  labs(x = "Vitamin Usage") +
  theme(text = element_text(size = 18))

pa2<-ggplot(Q3c.pred, aes(x = calories, y = alcohol)) +
  geom_point(size = 2)  +
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  theme(axis.title.y=element_blank())+
  labs(x = "Caloric Intake [Cal/day]") +
  theme(text = element_text(size = 18))

pa3<-ggplot(Q3c.pred, aes(x = fat, y = alcohol)) +
  geom_point(size = 2)  +
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  theme(axis.title.y=element_blank())+
  labs(x = "Fat Intake [mg/day]") +
  theme(text = element_text(size = 18))

pa4<-ggplot(Q3c.pred, aes(x = fiber, y = alcohol)) +
  geom_point(size = 2)  +
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  theme(axis.title.y=element_blank())+
  labs(x = "Dietary Fiber Intake [g/day]") +
  theme(text = element_text(size = 18))

pa5<-ggplot(Q3c.pred, aes(x = cholesterol, y = alcohol)) +
  geom_point(size = 2)  +
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  theme(axis.title.y=element_blank())+
  labs(x = "Cholesterol Intake [mg/day]") +
  theme(text = element_text(size = 18))


pa6<-ggplot(Q3c.pred, aes(x = betadiet, y = alcohol)) +
  geom_point(size = 2)  +
  geom_point(data = Q3c.pred[Alcoholic, ], 
             color = "blue", shape = 24, size = 3) +
  theme(axis.title.y=element_blank())+
  labs(x = "Dietary beta-carotene Intake [\U03BCg/day]") +
  theme(text = element_text(size = 18))


PlotsA<-ggarrange(pa1, pa2, pa3, pa4, pa5, pa6, labels="AUTO")
annotate_figure(PlotsA,
                top = text_grob("Alcoholic Intake vs other dietary variables", color = "black", face = "bold", size = 14),
                bottom = text_grob("The Blue triangle represents the extremely heavy drinker", color = "black",
                                   hjust = 1, x = 1, face = "italic", size = 10),
                left = text_grob("Alcohol Intake [L/week]", color = "black", rot = 90),
)
