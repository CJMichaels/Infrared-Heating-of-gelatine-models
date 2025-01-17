rm(list=ls())
setwd("~/Desktop/Reptile lamps publication/Final")

#Load necessary library (optional but recommended for handling CSV files)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyverse)
library(ggplot2)
library(survival)
library("survminer")

# Import the CSV file
data<-read.csv("dataoutlier.csv")

##1. Proportion of replicates that reached 25C for each basking lamp
#Pairwise Fisher's exact test
install.packages("fmsb")
library(fmsb)
successes <- c(39, 24, 37)
Total <- c(39, 38, 39)

pairwise.fisher.test(successes, Total)


##2. Linear models
#Choosing uncontrolled continuous explanatory variables 
#AIC testing of models with and without airtemp and air humidity 
AIC(lm(Duration~Lamp+Background+Distance+Humidity+Airtemp,data=data))#507.16
AIC(lm(Duration~Lamp+Background+Distance,data=data))#528.45

AIC(lm(Desiccation~Lamp+Background+Distance+Humidity+Airtemp,data=data)) #-1119.47
AIC(lm(Desiccation~Lamp+Background+Distance,data=data))#-1114.35
#Although models without air humidity and temperature had a lower AIC, biologically I find it necessary to include it in the linear model

AIC(lm(Ratio~Lamp+Background+Distance+Humidity+Airtemp,data=data))#114.46
AIC(lm(Ratio~Lamp+Background+Distance,data=data))#120.73

#Check for collinearity between air temperature and air humidity 
install.packages("usdm")
library(usdm)
myvars<-c ("Airtemp","Humidity")
vif(data[myvars]) 
cor.test(data$Airtemp,data$Humidity,method = "pearson")

#(a) Create Heating efficiency linear model
modelheat<-lm(Duration~Lamp+Background+Distance+Humidity+Airtemp,data=data)
summary(modelheat)

#Diagnostic plot
plot(modelheat)

#Remove outlier 
data2<-data[-c(35,37),]
modelheat<-lm(Duration~Lamp+Background+Distance+Humidity+Airtemp,data=data2)
summary(modelheat)
plot(modelheat)

# (b) Create surface to core temperature ratio model 
modelsc<-lm(Ratio~Lamp+Background+Distance+Humidity+Airtemp,data=data)
summary(modelsc)

#Diagnostic plot 
plot(modelsc)

#Remove outliers 
data3<-data[-c(37,91),]
modelsc<-lm(Ratio~Lamp+Background+Distance+Humidity+Airtemp,data=data3)
summary(modelsc)
plot(modelsc)

#(c) Create desiccation model 
modeld<-lm(Desiccation~Lamp+Background+Distance+Humidity+Airtemp,data=data)
summary(modeld)

#Diagnostic plot 
plot(modeld)

#Remover outlier 
data4<-data[-c(37),]
modeld<-lm(Desiccation~Lamp+Background+Distance+Humidity+Airtemp,data=data4)
summary(modeld)
plot(modeld)


##3.Pairwise post hoc tukey tests to look at how response variables differ among lamp types
install.packages("emmeans")
library(emmeans)
emmeans(modelheat, list(pairwise ~ Lamp), adjust = "tukey")
emmeans(modelsc, list(pairwise ~ Lamp), adjust = "tukey")
emmeans(modeld, list(pairwise ~ Lamp), adjust = "tukey")


##4.Calculate effect size of explanatory variables
install.packages("effectsize")
library(effectsize)
omega_squared(car::Anova(modelheat, type = 2),alternative = "two.sided")
omega_squared(car::Anova(modelsc, type = 2),alternative = "two.sided")
omega_squared(car::Anova(modeld, type = 2),alternative = "two.sided")


##5. Check whether air temperature differ among lamp type, background and Distance 
modelair<-lm(Airtemp~Lamp, data=data)
summary(modelair)
plot(modelair)
emmeans(modelair, list(pairwise ~ Lamp), adjust = "tukey")#confounding

modelairbackground<-lm(Airtemp~Background,data=data)
summary(modelairbackground)
plot(modelairbackground)#not confounding


modelairdistance<-lm(Airtemp~Distance, data=data)
summary(modelairdistance)
plot(modelairdistance)#confounding


##6. Check whether humidity differ among lamp type, background and Distance 
modelhumid<-lm(Humidity~Lamp,data=data)
summary(modelhumid)
plot(modelhumid)#not confounding

modelhumiditybackground<-lm(Humidity~Background, data=data)
summary(modelhumiditybackground)
plot(modelhumiditybackground)#not confounding 

modelhumiditydistance<-lm(Humidity~Distance, data=data)
summary(modelhumiditydistance)
plot(modelhumiditydistance)#not confounding


##7. Plots of raw data 
library(gridExtra)
ggplot(data, aes(x = Lamp, y = Duration,fill=Distance)) + labs(y="Heating efficiency (minutes)", x="Basking Lamp Types")+
  geom_boxplot(width=0.5)+theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(colour="black"),axis.text.y =  element_text(colour="black"),
        axis.title.x = element_text(vjust = -1.5),axis.title.y = element_text(vjust = +1.5))


plot1<- ggplot(data, aes(x = Lamp, y = Ratio,fill=Background)) + labs(y="Surface to core temperature change ratio", x="Basking Lamp Types")+
  geom_boxplot(width=0.4)+theme_bw()+theme(panel.border = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.background = element_blank(),
                                           axis.text.x = element_text(colour="black"),
                                           axis.text.y =  element_text(colour="black"),
                                           axis.line = element_line(color = "black"),
                                           axis.title.x = element_text(vjust = -1.5),
                                           axis.title.y = element_text(vjust = +1.5),
                                           plot.margin = unit (c(1,1,1,1),"cm")) 
#annotate("text", x=0.5,y=3,label = "(a)", hjust=0, vjust=-1.5,size=5, color = "black")

plot2<- ggplot(data, aes(x = Lamp, y = Desiccation)) + labs(y="Total dessication", x="Basking Lamp Types")+
  geom_boxplot(width=0.2,fill="salmon")+theme_bw()+theme(panel.border = element_blank(),
                                                         panel.grid.major = element_blank(),
                                                         panel.grid.minor = element_blank(),
                                                         panel.background = element_blank(),
                                                         axis.text.x = element_text(colour="black"),axis.text.y =  element_text(colour="black"),
                                                         axis.title.x = element_text(vjust = -1.5),axis.title.y = element_text(vjust = +1.5),
                                                         axis.line = element_line(color = "black"),
                                                         plot.margin = unit (c(1,1,1,1),"cm")) 
pdf("figs3.pdf",height=7,width=14)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()

