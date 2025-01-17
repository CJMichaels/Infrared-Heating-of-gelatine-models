
#Clear workspace
rm(list = ls())
setwd("~/Desktop/Reptile lamps publication/Chris R")

# Load necessary library (optional but recommended for handling CSV files)
library(readr)
library(tidyverse)
library(ggplot2)
library(survival)
library("survminer")


# Import the CSV file
sammi_data <- read_csv("sammi_data_with_dur.csv")


sammi_data_analysis <- sammi_data %>%
  select("Lamp", "Background", "Distance", "Replicate", "Duration", "Airtemp", "Humidity" = "Humidity(%)") %>%
  mutate(event = ifelse(Duration == 30, 0, 1)) %>%
  mutate(Airtemp = Airtemp- mean(Airtemp),
         Humidity = Humidity-mean(Humidity))

# Fit models

# Basic Cox PH model
cox_model_simple <- coxph(
  Surv((Duration), event) ~ Lamp + Background + Distance +  Airtemp + Humidity,
  data = sammi_data_analysis
)

# Test PH assumption
ph_test <- cox.zph(cox_model_simple)  # failed test; use AFT model as piece-wise and 
                                      # fractional polynomials likely inappropriate given size of dataset


# Fit AFT models with different error distributions
aft_model_normal <- survreg(Surv(Duration, event) ~ Distance + Airtemp + Humidity + Lamp + Background, 
                            data = sammi_data_analysis, dist = "gaussian")
aft_model_weibull <- survreg(Surv(Duration, event) ~ Distance + Airtemp + Humidity + Lamp + Background, 
                             data = sammi_data_analysis, dist = "weibull")
aft_model_lognormal <- survreg(Surv(Duration, event) ~ Distance + Airtemp + Humidity + Lamp + Background, 
                               data = sammi_data_analysis, dist = "lognormal")
aft_model_exponential <- survreg(Surv(Duration, event) ~ Distance + Airtemp + Humidity + Lamp + Background, 
                                 data = sammi_data_analysis, dist = "exponential")

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
plot(sammi_data_analysis$Duration, deviance_residuals, main = "Residuals vs Time", 
     xlab = "Duration", ylab = "Deviance Residuals", pch = 19)
abline(h = 0, col = "red")  # Zero line for reference

### Plots ###
# Get all unique combinations of Distance and Background
# Create an empty list to store KM plots
km_plots <- list()

# Get all unique combinations of Distance and Background
combinations <- unique(sammi_data_analysis[, c("Distance")])

# Loop through each combination
for (i in 1:nrow(combinations)) {
  subset_data <- subset(
    sammi_data_analysis,
    Distance == combinations$Distance[i] 
  )
  # factor() make sure "halogen" comes before "cermaics"
  subset_data$Lamp <- factor (subset_data$Lamp, levels = c("Halogen",
                                                          "Deep Heat Projector",
                                                           "Ceramic"))
  # Check if there's enough data for the subset
  if (nrow(subset_data) > 0) {
    # Fit KM survival curves for the subset
    km_fit <- survfit(Surv(Duration, event) ~ Lamp, data = subset_data)
    strata_labels <- gsub("Lamp=", "", names(km_fit$strata))
    strata_labels <- gsub("Deep Heat Projector", "Carbon Filament Heater", strata_labels)
    
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

km_plots
class(km_plots) #list
plot2a <- km_plots[[1]]
plot2b <- km_plots[[2]]
plot2a <-recordPlot({plot2a})
plot2b <-recordPlot({plot2b})
my_plots <- list (plot2a, plot2b)
library(gridExtra)
pdf("figs2.pdf",height=8,width=10)
for (i in seq_along(km_plots)){
  print(km_plots[[i]])
}
dev.off()

# when u increase the pixel (100,300,600), the actual image dimensions in pixels 
#will also increase if you want the physical size (in inches) to remain constant
#--> need to multiply the dimensions by the pixcel
for (i in 1:2){
  jpeg(paste0("figs2_",i,".jpeg"), width = 11*200, height = 10*200, units = "px", res=200)
  print(km_plots[[1]])
  print(km_plots[[2]])
}
dev.off()

