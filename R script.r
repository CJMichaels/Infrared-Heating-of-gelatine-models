rm(list=ls())

#Load necessary library (optional but recommended for handling CSV files)
library(readr)
library(dplyr)
library(ggplot2)
library(survival)
library("survminer")

# Import the CSV file
data<-read.csv("data_R.csv")

# Set the reference level as Ceramic lamp
data$Lamp <- as.factor(data$Lamp)
data$Lamp_type <- relevel(data$Lamp, ref = "Ceramic")

################################################################################
##1. Proportion of replicates that reached 25C for each basking lamp
#Pairwise Fisher's exact test
library(fmsb)
successes <- c(39, 24, 37)
Total <- c(39, 38, 39)

pairwise.fisher.test(successes, Total)

################################################################################
##2. Survival models for heating effectiveness
data_analysis <- data[, c("Lamp_type", "Background", "Distance", "Replicate", "Duration", "Airtemp", "Humidity")] 
data_analysis <- data_analysis %>%
  mutate(event = ifelse(Duration == 30, 0, 1)) %>%
  mutate(Airtemp = Airtemp- mean(Airtemp),
         Humidity = Humidity-mean(Humidity))

data_analysis$Duration <- as.numeric(data_analysis$Duration)

#Fit models
#Basic Cox PH model
cox_model_simple <- coxph(
  Surv((Duration), event) ~ Lamp_type + Background + Distance +  Airtemp + Humidity,
  data = data_analysis)

# Test PH assumption
ph_test <- cox.zph(cox_model_simple) #failed test; use AFT model as piece-wise and 
                                     #fractional polynomials likely inappropriate given size of dataset


# Fit AFT models with different error distributions
aft_model_normal <- survreg(Surv(Duration, event) ~ Distance + Airtemp + Humidity + Lamp_type + Background, 
                            data = data_analysis, dist = "gaussian")
aft_model_weibull <- survreg(Surv(Duration, event) ~ Distance + Airtemp + Humidity + Lamp_type + Background, 
                             data = data_analysis, dist = "weibull")
aft_model_lognormal <- survreg(Surv(Duration, event) ~ Distance + Airtemp + Humidity + Lamp_type + Background, 
                               data = data_analysis, dist = "lognormal")
aft_model_exponential <- survreg(Surv(Duration, event) ~ Distance + Airtemp + Humidity + Lamp_type + Background, 
                                 data = data_analysis, dist = "exponential")


# Compare models using AIC
aic_values <- AIC(aft_model_normal, aft_model_weibull, aft_model_lognormal, aft_model_exponential)

# Compare models using BIC
bic_values <- BIC(aft_model_normal, aft_model_weibull, aft_model_lognormal, aft_model_exponential)


# Lognormal best fit
# Check model summaries (coefficients, standard errors, and significance)
summary(aft_model_lognormal)

# Deviance residuals
deviance_residuals <- residuals(aft_model_lognormal, type = "deviance")
hist(deviance_residuals, main = "Histogram of Deviance Residuals", xlab = "Deviance Residuals", col = "lightblue", breaks = 20)

# dfbetas (influence of each predictor on the coefficient)
dfbetas_values <- residuals(aft_model_lognormal, type = "dfbetas")
hist(dfbetas_values)

# Check residuals for normality
qqnorm(deviance_residuals)
qqline(deviance_residuals, col = "red")

# Residuals vs Time plot
plot(data_analysis$Duration, deviance_residuals, main = "Residuals vs Time", 
     xlab = "Duration", ylab = "Deviance Residuals", pch = 19)
abline(h = 0, col = "red")  # Zero line for reference


# Plots for survival model
# Get all unique combinations of Distance and Background
# Create an empty list to store KM plots
km_plots <- list()

# Get all unique combinations of Distance and Background
combinations <- data_analysis %>% distinct(Distance)
#combinations <- data.frame(Distance = c("20cm", "30cm"))

# Loop through each combination
for (i in 1:nrow(combinations)) {
  subset_data <- subset(
    data_analysis,
    Distance == combinations$Distance[i] 
  )
  
# factor() make sure "Incandescent" comes before "ceramics"
subset_data$Lamp <- factor (subset_data$Lamp, 
                               levels = c("Incandescent", "Carbon Filament Heater", "Ceramic"))
# Check if there's enough data for the subset
  if (nrow(subset_data) > 0) {
    # Fit KM survival curves for the subset
    km_fit <- survfit(Surv(Duration, event) ~ Lamp, data = subset_data)
    strata_labels <- gsub("Lamp=", "", names(km_fit$strata))
    
    # Generate KM plot for the subset
    plot_title <- paste(
      "Distance:", combinations$Distance[i]
    )
    plot <- ggsurvplot(
      km_fit,
      data = subset_data,
      xlab = "Time (minutes)",
      ylab = "Probability of not being at core target temperature",
      title = plot_title,
      legend.title = "",
      legend.labs = strata_labels,
      palette = "Dark2",
      risk.table = TRUE,                  # Show the numbers at risk table
      risk.table.col = "strata",           # Optional: color the risk table by strata (group)
      risk.table.height = 0.25,             # Adjust the height of the risk table
      risk.table.y.text.col = TRUE,        # Color the numbers at risk
      risk.table.fontsize = 3#,             # Adjust font size of numbers in the table
      # risk.table.title = "Number at Risk"  # Title of the numbers at risk table
    )
    
# Save the plot to the list
  km_plots[[plot_title]] <- plot
  } else {
    message(paste(
      "Not enough data for Distance =", combinations$Distance[i], 
      "and Background =", combinations$Background[i]
    ))
  }
}

km_plots #list

#Export the plot in jpeg format
for (i in 1:2){
  jpeg(paste0("figs2a_",i,".jpeg"), width = 11*200, height = 10*200, units = "px", res=200)
  print(km_plots[[1]])
  print(km_plots[[2]])
}
dev.off()
################################################################################
##3. Linear models for desiccation and surface to core temperature ratio

# Choosing uncontrolled continuous explanatory variables 
# AIC testing of models with and without airtemp and air humidity 
AIC(lm(Desiccation~Lamp+Background+Distance+Humidity+Airtemp,data=data)) #-1119.47
AIC(lm(Desiccation~Lamp+Background+Distance,data=data))#-1114.35
#Although models without air humidity and temperature had a lower AIC, biologically I find it necessary to include it in the linear model

AIC(lm(Ratio~Lamp+Background+Distance+Humidity+Airtemp,data=data))#114.46
AIC(lm(Ratio~Lamp+Background+Distance,data=data))#120.73

# Check for collinearity between air temperature and air humidity 
library(usdm)
myvars<-c ("Airtemp","Humidity")
vif(data[myvars]) 
cor.test(data$Airtemp,data$Humidity,method = "pearson")


# (a)Create surface to core temperature ratio model 
modelsc<-lm(Ratio~Lamp_type+Background+Distance+Humidity+Airtemp,data=data)
summary(modelsc)

# Diagnostic plot 
plot(modelsc)

# Remove outliers 
data3<-data[-c(37,91),]
modelsc<-lm(Ratio~Lamp_type+Background+Distance+Humidity+Airtemp,data=data3)
summary(modelsc)
plot(modelsc)


#(b) Create desiccation model 
modeld<-lm(Desiccation~Lamp_type+Background+Distance+Humidity+Airtemp,data=data)
summary(modeld)

# Diagnostic plot 
plot(modeld)

# Remover outlier 
data4<-data[-c(37),]
modeld<-lm(Desiccation~Lamp_type+Background+Distance+Humidity+Airtemp,data=data4)
summary(modeld)
plot(modeld)


# Pairwise post hoc tukey tests to look at how response variables differ among lamp types
library(emmeans)
emmeans(modelsc, list(pairwise ~ Lamp_type), adjust = "tukey")
emmeans(modeld, list(pairwise ~ Lamp_type), adjust = "tukey")


# Calculate effect size of explanatory variables
library(effectsize)
omega_squared(car::Anova(modelsc, type = 2),alternative = "two.sided")
omega_squared(car::Anova(modeld, type = 2),alternative = "two.sided")


# Check whether air temperature differ among lamp type, background and Distance 
modelair<-lm(Airtemp~Lamp_type, data=data)
summary(modelair)
plot(modelair)
emmeans(modelair, list(pairwise ~ Lamp_type), adjust = "tukey")#confounding

modelairbackground<-lm(Airtemp~Background,data=data)
summary(modelairbackground)
plot(modelairbackground)#not confounding

modelairdistance<-lm(Airtemp~Distance, data=data)
summary(modelairdistance)
plot(modelairdistance)#confounding


# Check whether humidity differ among lamp type, background and Distance 
modelhumid<-lm(Humidity~Lamp,data=data)
summary(modelhumid)
plot(modelhumid)#not confounding

modelhumiditybackground<-lm(Humidity~Background, data=data)
summary(modelhumiditybackground)
plot(modelhumiditybackground)#not confounding 

modelhumiditydistance<-lm(Humidity~Distance, data=data)
summary(modelhumiditydistance)
plot(modelhumiditydistance)#not confounding


# Plots of linear model
library(gridExtra)
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

plot2<- ggplot(data, aes(x = Lamp, y = Desiccation)) + labs(y="Total dessication", x="Basking Lamp Types")+
  geom_boxplot(width=0.2,fill="salmon")+theme_bw()+theme(panel.border = element_blank(),
                                                         panel.grid.major = element_blank(),
                                                         panel.grid.minor = element_blank(),
                                                         panel.background = element_blank(),
                                                         axis.text.x = element_text(colour="black"),axis.text.y =  element_text(colour="black"),
                                                         axis.title.x = element_text(vjust = -1.5),axis.title.y = element_text(vjust = +1.5),
                                                         axis.line = element_line(color = "black"),
                                                       plot.margin = unit (c(1,1,1,1),"cm")) 


# Export the plot in jpeg format
library(cowplot)
combined_plot <- plot_grid(plot1, plot2, ncol=2)
ggsave("fig3.jpeg", combined_plot, width = 13, height = 5)
