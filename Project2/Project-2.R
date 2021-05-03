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
plasma$plasmacat <- factor(plasma$lowplasma,
                           levels = c(0, 1),
                           labels = c("high", "low"))
# How many are low?
low_total <- sum(PositivePlasma$lowplasma)
high_total <- 315 - low_total
low_percentage <- (low_total/316)*100 

table(PositivePlasma$smokstat, PositivePlasma$lowplasma)



































