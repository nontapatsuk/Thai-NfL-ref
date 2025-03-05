## *******************************************************************************************************************************************
## Disclaimer: This code was developed for research purposes and may not be optimised for clinical or commercial applications. 
## Users should validate the results independently before relying on them.
## *******************************************************************************************************************************************

## Project: Normal sNfL in Thai adults
## Purpose: This script runs linear regression model and GAMLSS
## Author: Nontapat Sukhonpanich, MD                                              
## nontapat.suk@mahidol.ac.th                                        
## March 3rd, 2025
## Division of Neurology, Department of Medicine, Faculty of Medicine Siriraj Hospital, Mahidol University
## Bangkok, Thailand
## *******************************************************************************************************************************************
###################################################################################################################
## ---------------------------------------------------------------------------------------------------------------

rm(list = ls())
closeAllConnections()
options(stringsAsFactors = FALSE)

library(tidyverse)
library(data.table)
library(xlsx)
library(nlme)
library(ggplot2)
library(dplyr)
library(ggpubr)
# install.packages("gamlss")
library(gamlss)
# install.packages("boot")
library(boot)
options(scipen=999)

library(readxl)

#Load dataset
setwd("D:/Si_Neuro/NfL normal") ##choose your own directory
df <- read_excel("Normal sNfL proj_BB_22.8.24.xlsx", 
                 col_types = c("text", "numeric", "text", 
                               "text", "text", "text", "numeric"))
View(df) 
str(df)

#Change to factor
df$ID <- as.factor(df$ID)
df$Sex_1F0M <- as.factor(df$Sex_1F0M)
df$Age_Cat <- as.factor (df$Age_Cat) #1=20-29, 2=30-29, 3=40-49, 4=50-59, 5=60-69
levels(df$Age_Cat)<-c("20-29","30-39","40-49","50-59","60-69")

#Check distribution of sNfL
hist(df$sNfL) #left skew
qqnorm (df$sNfL)
qqline (df$sNfL, col = "red")

shapiro_test_result <- shapiro.test(df$sNfL)
print(shapiro_test_result) #p=0.0000000001894

#clear NA
df_a <- subset(df, Age!=is.na(Age))
df2 <- subset(df_a, Sex_1F0M!=is.na(Sex_1F0M))
df2 <- subset(df2, sNfL!=0.53)
View(df2) #n=223

hist (df2$Age)
shapiro_test_result2 <- shapiro.test(df2$Age)

#Serum NfL - summarised
df2 %>% 
  summarise(n(), mean = mean (sNfL),sd=sd(sNfL), median=median(sNfL), 
            P2.5=quantile(sNfL,0.025),P97.5=quantile(sNfL,0.975),
            min=min(sNfL),max=max(sNfL))
#`n()`  mean    sd median  P2.5 P97.5   min   max
#<int> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
# 223  6.57  3.36   5.83  2.29  15.9     1  18.4

boxplot(df2$sNfL, horizontal = TRUE,
        xlab = "Serum NfL (pg/mL)",
        ylab = "Age(years)")

help ("boxplot")

##median and 95%CI
df_median_all <- df2$sNfL

median_func <- function(df_median_all, indices) {
  sampled_data <- df_median_all[indices]
  return(median(sampled_data))
}

bootstrap_results_all <- boot(df_median_all, statistic = median_func, R = 1000) # R is the number of resamples

ci_all <- boot.ci(bootstrap_results, type = "perc")
print(ci_all)

#95% CI of median by sex
df_male <- df2 %>% 
  filter(Sex_1F0M==0)
df_male_nfl <- df_male$sNfL

median_func_m <- function(df_male_nfl, indices) {
  sampled_data <- df_male_nfl[indices]
  return(median(sampled_data))
}

bootstrap_results_m <- boot(df_male_nfl, statistic = median_func_m, R = 1000) # R is the number of resamples

ci_m <- boot.ci(bootstrap_results_m, type = "perc")
print(ci_m)

#-#
df_female <- df2 %>% 
  filter(Sex_1F0M==1)
df_female_nfl <- df_female$sNfL

median_func_fm <- function(df_female_nfl, indices) {
  sampled_data <- df_female_nfl[indices]
  return(median(sampled_data))
}

bootstrap_results_fm <- boot(df_female_nfl, statistic = median_func_fm, R = 1000) # R is the number of resamples

ci_fm <- boot.ci(bootstrap_results_fm, type = "perc")
print(ci_fm)

## difference between male and female
wilcox.test(df_male_nfl, df_female_nfl)

##by age category
#20-29
df_20 <- df2 %>% 
  filter(Age_Cat=="20-29")
df_20_nfl<- df_20$sNfL

median_func_20 <- function(df_20_nfl, indices) {
  sampled_data <- df_20_nfl[indices]
  return(median(sampled_data))
}

bootstrap_results_20 <- boot(df_20_nfl, statistic = median_func_20, R = 1000) # R is the number of resamples

ci_20 <- boot.ci(bootstrap_results_20, type = "perc")
print(ci_20)

#30-39
df_30 <- df2 %>% 
  filter(Age_Cat=="30-39")
df_30_nfl<- df_30$sNfL

median_func_30 <- function(df_30_nfl, indices) {
  sampled_data <- df_30_nfl[indices]
  return(median(sampled_data))
}

bootstrap_results_30 <- boot(df_30_nfl, statistic = median_func_30, R = 1000) # R is the number of resamples

ci_30 <- boot.ci(bootstrap_results_30, type = "perc")
print(ci_30)

#40-49
df_40 <- df2 %>% 
  filter(Age_Cat=="40-49")
df_40_nfl<- df_40$sNfL

median_func_40 <- function(df_40_nfl, indices) {
  sampled_data <- df_40_nfl[indices]
  return(median(sampled_data))
}

bootstrap_results_40 <- boot(df_40_nfl, statistic = median_func_40, R = 1000) # R is the number of resamples

ci_40 <- boot.ci(bootstrap_results_40, type = "perc")
print(ci_40)

#50-59
df_50 <- df2 %>% 
  filter(Age_Cat=="50-59")
df_50_nfl<- df_50$sNfL

median_func_50 <- function(df_50_nfl, indices) {
  sampled_data <- df_50_nfl[indices]
  return(median(sampled_data))
}

bootstrap_results_50 <- boot(df_50_nfl, statistic = median_func_50, R = 1000) # R is the number of resamples

ci_50 <- boot.ci(bootstrap_results_50, type = "perc")
print(ci_50)

#60-69
df_60 <- df2 %>% 
  filter(Age_Cat=="60-69")
df_60_nfl<- df_60$sNfL

median_func_60 <- function(df_60_nfl, indices) {
  sampled_data <- df_60_nfl[indices]
  return(median(sampled_data))
}

bootstrap_results_60 <- boot(df_60_nfl, statistic = median_func_60, R = 1000) # R is the number of resamples

ci_60 <- boot.ci(bootstrap_results_60, type = "perc")
print(ci_60)

#20-39
df2$Age_Cat2 <- if_else(df2$Age_Cat=="20-29"|df2$Age_Cat=="30-39","20-39","40-69")

df_2039 <- df2 %>% 
  filter(Age_Cat2=="20-39")
df_2039_nfl<- df_2039$sNfL

median_func_2039 <- function(df_2039_nfl, indices) {
  sampled_data <- df_2039_nfl[indices]
  return(median(sampled_data))
}

bootstrap_results_2039 <- boot(df_2039_nfl, statistic = median_func_2039, R = 1000) # R is the number of resamples

ci_2039 <- boot.ci(bootstrap_results_2039, type = "perc")
print(ci_2039)

#40-69
df_4069 <- df2 %>% 
  filter(Age_Cat2=="40-69")
df_4069_nfl<- df_4069$sNfL

median_func_4069 <- function(df_4069_nfl, indices) {
  sampled_data <- df_4069_nfl[indices]
  return(median(sampled_data))
}

bootstrap_results_4069 <- boot(df_4069_nfl, statistic = median_func_4069, R = 1000) # R is the number of resamples

ci_4069 <- boot.ci(bootstrap_results_4069, type = "perc")
print(ci_4069)

#explore each cat
sNfL_age_cat_nona <- df2 %>% 
  group_by(Age_Cat) %>% 
  summarise(n(), mean = mean(sNfL), sd = sd(sNfL), median=median(sNfL), 
            P2.5=quantile(sNfL,0.025),P97.5=quantile(sNfL,0.975),
            min=min(sNfL),max=max(sNfL))
sNfL_age_cat_nona

bp_sNfL_age_cat <- boxplot(df2$sNfL ~ df2$Age_Cat, horizontal = TRUE,
                           xlab = "Serum NfL (pg/mL)",
                           ylab = "Age Category")

df2$Age_Cat2 <- if_else(df2$Age_Cat=="20-29"|df2$Age_Cat=="30-39","20-39","40-69")

sNfL_age_cat2_nona <- df2 %>% 
  group_by(Age_Cat2) %>% 
  summarise(n(), mean = mean(sNfL), sd= sd(sNfL), 
            median=median(sNfL), P2.5=quantile(sNfL,0.025),
            P97.5=quantile(sNfL,0.975),min=min(sNfL),max=max(sNfL))
sNfL_age_cat2_nona

bp_sNfL_age_cat2_nona <- boxplot(df2$sNfL ~ df2$Age_Cat2, horizontal = TRUE,
                            xlab = "Serum NfL (pg/mL)",
                            ylab = "Age Category")

#group by sex
sNfL_sex_nona <- df2 %>% 
  group_by(Sex_1F0M) %>% 
  summarise(n(), mean = mean(sNfL), sd= sd(sNfL), 
            median=median(sNfL), P2.5=quantile(sNfL,0.025),
            P97.5=quantile(sNfL,0.975),min=min(sNfL),max=max(sNfL))
sNfL_sex_nona

bp_sNfL_sex_nona <- boxplot(df2$sNfL ~ df2$Sex_1F0M, horizontal = TRUE,
                                 xlab = "Serum NfL (pg/mL)",
                                 ylab = "Sex (0M/1F)")

#####--------------------------------------------------------------#####

##Correlation between log NfL and age + sex
#log transform
df2$logNfL <- log10(df2$sNfL)

#create an upper 95% prediction line
model_log1 <- lm(logNfL ~ Age, data=df2)
summary (model_log1) #adjusted R square 0.405

###Equation for only age ( no sex adjustment )
### log sNfL = 0.27237 + 0.01239(Age)
residual_age <- residuals(model_log1)
qqnorm (residual_age) #normal QQ plot
hist(residual_age) #normal histogram

summary(residual_age)
sd (residual_age) #0.1774778

###Upper limit (mean + 1.96 SD)
(1.96*0.1774778)
###log sNfL = 0.27237 + 0.01239(Age) + 0.3478565

df_sNfL_age <- data.frame(Age = c(20, 30, 40, 50, 60, 70))
df_sNfL_age$Cal_sNfL <- 10^(0.27237 + 0.01239*df_sNfL_age$Age + 0.3478565)
view(df_sNfL_age)

###Log10 sNfL  and sex
model2_sex <- lm (logNfL ~ Sex_1F0M, data=df2)
summary(model2_sex) #p=0.0206, c = 0.80015, coef = -0.07191

##Full model - sNfL by age and adjusted for sex
model_log <- lm(logNfL ~ Age + Sex_1F0M, data=df2)
summary (model_log) 
print(model_log)
residual2 <- residuals(model_log)
sd (residual2) #0.174377

###equation for predicted mean: log(10) sNfL = (0.012339)x(Age)+(-0.066608 )x(Sex)+0.312393
df2$Sex_num = as.numeric(df2$Sex_1F0M)

df_sNfL_age$Cal_sNfL_male <- 10^(0.012006 * df_sNfL_age$Age + (-0.059990) * 0 + 0.325539 + 1.96*0.1657)
df_sNfL_age$Cal_sNfL_female <- 10^(0.012006 * df_sNfL_age$Age + (-0.059990) * 1 + 0.325539  + 1.96*0.1657)

View(df_sNfL_age)

# Create a scatter plot with regression line
q <- (ggscatter(
  df2,
  x = "Age",
  y = "logNfL",
  color = "Sex",  # Color by Sex
  shape = "Sex",        
  add = "reg.line",
  conf.int = FALSE,
  size = 1.4,
  cor.coef = FALSE,
  cor.method = "pearson",
  xlab = "Age (year)",
  ylab = "Log sNfL",
  ylim = c(0, 1.5),
  add.params = list(size = 0.7)
)) +
  scale_color_manual(
    name = "Sex",             # Change legend title
    labels = c("Female", "Male"), # Change legend labels
    values = c("#F05039", "#3D65A5") # Define custom colors
  ) +
  scale_shape_manual(
    name = "Sex",
    labels = c("Female", "Male"),
    values = c(16, 17)  # Choose different shapes (16 = filled circle, 17 = filled triangle)
  ) +
  stat_cor(
    aes(label = ..rr.label..),
    label.x.npc = "left",
    label.y.npc = "top",
    output.type = "text",
    method = "pearson"
  )

# Display the plot
print(q)

# Save the plot --> figure 1
ggexport(q, width=4000 ,height=2500, filename = "LM_NfL_age_sex.tiff", res = 550)

# Plot for graphical abstract
q2 <- (ggscatter(
  df2,
  x = "Age",
  y = "logNfL",
  color = "Sex",  # Color by Sex
  shape = "Sex",        
  add = "reg.line",
  conf.int = FALSE,
  size = 2.0,
  cor.coef = FALSE,
  cor.method = "pearson",
  xlab = "Age (year)",
  ylab = "Log sNfL",
  ylim = c(0, 1.5),
  add.params = list(size = 0.8)
)) +
  scale_color_manual(
    name = "Sex",             # Change legend title
    labels = c("Female", "Male"), # Change legend labels
    values = c("#F05039", "#3D65A5") # Define custom colors
  ) +
  scale_shape_manual(
    name = "Sex",
    labels = c("Female", "Male"),
    values = c(16, 17)  # Choose different shapes (16 = filled circle, 17 = filled triangle)
  ) +
  stat_cor(
    aes(label = ..rr.label..),
    label.x.npc = "left",
    label.y.npc = "top",
    output.type = "text",
    method = "pearson"
  )

ggexport(q2, width=3000 ,height=2100, filename = "LM_NfL_age_sex2.tiff", res = 550)

#####-----------------------------------------------------------------#####
###Modeling using GAMLSS
###z-score
install.packages("gamlss") #if not done
library(gamlss)

df3 <- df2 %>% 
  select(sNfL,Age,Sex_1F0M)
View (df3)

# Fit a GAMLSS model to the data
## model for age
model <- gamlss(sNfL ~ pb(Age),       
                sigma.formula = ~ pb(Age), 
                family = BCCG,            
                data = df3)

summary(model)

## model for age and sex
model2 <- gamlss(sNfL ~ pb(Age)+Sex_1F0M,       # Smoothing function for age
                sigma.formula = ~ pb(Age), # Allow scale to vary with age
                family = BCCG,            # Choose an appropriate family
                data = df3)

# Diagnostics
wp (model2)
plot (model2)
summary (model2)

df5 <- df3
df5$muNfL_s <- fitted(model2, "mu")
df5$sigmaNfL_s <- fitted(model2, "sigma")
df5$nuNfL_s <- fitted(model2,"nu")

df5$z_score <- ((df5$sNfL /df5$muNfL_s )^df5$nuNfL_s - 1) / (df5$sigmaNfL_s  * df5$nuNfL_s)
hist(df5$z_score)
sum(df5$z_score > 1.96 | df5$z_score < -1.96, na.rm = TRUE)
sum(df5$z_score > 1.96, na.rm = TRUE)

#### Z-score
plot_data_sex <- data.frame(Age = df5$Age, Z_Score = df5$z_score)

Z_sex <- ggplot(plot_data_sex, aes(x = Age, y = Z_Score, 
                                   color = df5$Sex_1F0M, 
                                   shape = df5$Sex_1F0M)) +  # Map shape to Sex
  geom_point() +
  geom_hline(yintercept = c(-1.96, 0, 1.96), linetype = "dashed", color = "red") +
  scale_color_manual(name = "Sex",             
                     labels = c("Female", "Male"),
                     values = c("0" = "#F05039", "1" = "#3D65A5")) +  # Custom colors
  scale_shape_manual(name = "Sex",
                     labels = c("Female", "Male"),
                     values = c(16, 17)) +  # 16 = Circle, 17 = Triangle
  labs(title = "Age-Normative Z-Scores",
       x = "Age (year)",
       y = "Z-Score") +
  theme_minimal()

# save the plot --> supplementary fig 1
ggexport(Z_sex, width=4500 ,height=3000, filename = "Z_sex.tiff", res = 550)

## Centiles
# Define the 97.5th percentile quantile
q97_5 <- 0.975

# Create a sequence of ages for prediction
age_seq <- seq(min(df3$Age), max(df3$Age), length.out = 100)

# Compute the 97.5th percentile (total)
# Predict μ, σ, and ν at
mu <- predict(model, newdata = data.frame(Age = age_seq), what = "mu", type = "response")
sigma <- predict(model, newdata = data.frame(Age = age_seq), what = "sigma", type = "response")
nu <- predict(model, newdata = data.frame(Age = age_seq), what = "nu", type = "response")

P97.5 <- qBCCG(0.975, mu = mu, sigma = sigma, nu = nu)

centiles(model, xvar = df3$Age, cent = c(2.5, 10, 50, 90, 97.5))
centiles(model, xvar = df3$Age, cent = 97.5)

plot_data <- data.frame(Age = age_seq, P97.5 = P97.5)

# Create a prediction data frame for both sexes (97.5th percentile)
# Create a grid of Age values for prediction
age_grid <- seq(min(df3$Age), max(df3$Age), length.out = 100)

predict_data3 <- data.frame(
  Age = rep(age_grid, 2),
  Sex_1F0M = c(rep(0, length(age_grid)), rep(1, length(age_grid)))  # 0 = Male, 1 = Female
)

mu_sex2 <- predict(model2, newdata = predict_data3, what = "mu", type = "response")
sigma_sex2 <- predict(model2, newdata = predict_data3, what = "sigma", type = "response")
nu_sex2 <- predict(model2, newdata = predict_data3, what = "nu", type = "response")

predict_data3$P97.5 <- qBCCG(0.975, mu = mu_sex2, sigma = sigma_sex2, nu = nu_sex2)

# Plot the percentiles
library(ggplot2)

P97.5_by_age_sex <- ggplot(df3, aes(x = Age, y = sNfL)) +
  geom_point(aes(color = factor(Sex_1F0M),shape=factor(Sex_1F0M))) +  # Raw data points
  geom_line(data = plot_data, aes(x = Age, y = P97.5), 
            color = "red", size = 0.7) +    
  geom_line(data = predict_data3, 
            aes(x = Age, y = P97.5, color = factor(Sex_1F0M)), 
            size = 0.7, linetype="dashed") +                                    # 97.5th percentile lines
  labs(title = "Age-Normative 97.5th Percentile Curves for Serum NfL by Sex",
       x = "Age (years)",
       y = "Serum NfL Level (pg/mL)",
       color = "Sex") +                                   # Change legend title
  scale_shape_manual(name = "Sex",
                     labels = c("Male", "Female"),
                     values = c(17, 16)) +  # 16 = Circle, 17 = Triangle
  scale_color_manual(                                     # Change legend categories
    values = c("0" = "#3D65A5", "1" = "#F05039"),              # Set custom colors
    labels = c("Male", "Female")                         # Change category labels
  ) +
  theme_minimal()

#save --> figure 3
ggexport(P97.5_by_age_sex, width=4500 ,height=3000, filename = "sNfL_P97.5_by_age_sex.tiff", res = 550)

#table supplementary 1
age_new <- c(20,25,30,35,40,45,50,55,60,65)
predict_data4 <- data.frame(
  Age = rep(age_new, 2),
  Sex_1F0M = c(rep(0, length(age_new)), rep(1, length(age_new)))  # 0 = Male, 1 = Female
)

mu_new2 <- predict(model2, newdata = predict_data4, what = "mu", type = "response")
sigma_new2<- predict(model2, newdata = predict_data4, what = "sigma", type = "response")
nu_new2 <- predict(model2, newdata = predict_data4, what = "nu", type = "response")

predict_data4$P97.5 <- qBCCG(0.975, mu = mu_new2, sigma = sigma_new2, nu = nu_new2)
predict_data4$P2.5 <- qBCCG(0.025, mu = mu_new2, sigma = sigma_new2, nu = nu_new2)
predict_data4$P5 <- qBCCG(0.05, mu = mu_new2, sigma = sigma_new2, nu = nu_new2)
predict_data4$P10 <-qBCCG(0.1, mu = mu_new2, sigma = sigma_new2, nu = nu_new2)
predict_data4$P50 <- qBCCG(0.5, mu = mu_new2, sigma = sigma_new2, nu = nu_new2)
predict_data4$P90 <- qBCCG(0.9, mu = mu_new2, sigma = sigma_new2, nu = nu_new2)
predict_data4$P95 <- qBCCG(0.95, mu = mu_new2, sigma = sigma_new2, nu = nu_new2)
View(predict_data4)

#Centile 2.5,5,10,50,90,95,97.5
P97.5 <- qBCCG(0.975, mu = mu, sigma = sigma, nu = nu)
P95 <- qBCCG(0.95, mu = mu, sigma = sigma, nu = nu)
P90 <- qBCCG(0.90, mu = mu, sigma = sigma, nu = nu)
P50 <- qBCCG(0.50, mu = mu, sigma = sigma, nu = nu)
P10 <- qBCCG(0.10, mu = mu, sigma = sigma, nu = nu)
P5 <- qBCCG(0.05, mu = mu, sigma = sigma, nu = nu)
P2.5 <- qBCCG(0.025, mu = mu, sigma = sigma, nu = nu)

plot_data_centile <- data.frame(Age = age_seq,P2.5 = P2.5, P5=P5, P10=P10,
                                P50=P50, P90=P90, P95=P95,P97.5 = P97.5)


#new fig centiles
P_by_age <- ggplot(df3, aes(x = Age, y = sNfL)) +
  geom_point(aes(color = factor(Sex_1F0M),shape=factor(Sex_1F0M))) +  # Raw data points
  geom_line(data = plot_data_centile, aes(x = Age, y = P97.5), 
            color = "#f94144", linewidth = 0.5) +  
  geom_line(data = plot_data_centile, aes(x = Age, y = P2.5), 
            color = "#577590", linewidth = 0.5) +  
  geom_line(data = plot_data_centile, aes(x = Age, y = P5), 
            color = "#43aa8b", linewidth = 0.5) +  
  geom_line(data = plot_data_centile, aes(x = Age, y = P10), 
            color = "#90be6d", linewidth = 0.5) + 
  geom_line(data = plot_data_centile, aes(x = Age, y = P50), 
            color = "#f9c74f", linewidth = 0.5) +  
  geom_line(data = plot_data_centile, aes(x = Age, y = P90), 
            color = "#f8961e", linewidth = 0.5) +  
  geom_line(data = plot_data_centile, aes(x = Age, y = P95), 
            color = "#f3722c", linewidth = 0.5) +  
  labs(title = "Age-Normative Percentile Curves for Serum NfL",
       x = "Age (years)",
       y = "Serum NfL Level (pg/mL)",
       color = "Sex") +                                   # Change legend title
  scale_shape_manual(name = "Sex",
                     labels = c("Male", "Female"),
                     values = c(17, 16)) +  # 16 = Circle, 17 = Triangle
  scale_color_manual(                                     # Change legend categories
    values = c("0" = "#3D65A5", "1" = "#F05039"),              # Set custom colors
    labels = c("Male", "Female")                         # Change category labels
  ) +
  # Add labels for each curve
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)),  # Add labels at the max Age
            aes(x = Age, y = P97.5, label = "97.5"), 
            hjust = -0.2, color = "#f94144",size=3) +
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)), 
            aes(x = Age, y = P2.5, label = "2.5"), 
            hjust = -0.2, color = "#577590",size=3) +
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)), 
            aes(x = Age, y = P5, label = "5"), 
            hjust = -0.3, color = "#43aa8b",size=3) +
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)), 
            aes(x = Age, y = P10, label = "10"), 
            hjust = -0.2, color = "#90be6d",size=3) +
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)), 
            aes(x = Age, y = P50, label = "50"), 
            hjust = -0.2, color = "#f9c74f",size=3) +
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)), 
            aes(x = Age, y = P90, label = "90"), 
            hjust = -0.3, color = "#f8961e",size=3) +
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)), 
            aes(x = Age, y = P95, label = "95"), 
            hjust = -0.2, color = "#f3722c",size=3) +
  theme_minimal()

#save --> figure 2
ggexport(P_by_age, width=4500 ,height=3000, filename = "percentile_curves.tiff", res = 550)

# Graphical abstract fig 2
P_by_age2 <- ggplot(df3, aes(x = Age, y = sNfL)) +
  geom_point(aes(color = factor(Sex_1F0M),shape=factor(Sex_1F0M)),size=2) +  # Raw data points
  geom_line(data = plot_data_centile, aes(x = Age, y = P97.5), 
            color = "#f94144", linewidth = 0.5) +  
  geom_line(data = plot_data_centile, aes(x = Age, y = P2.5), 
            color = "#577590", linewidth = 0.5) +  
  geom_line(data = plot_data_centile, aes(x = Age, y = P5), 
            color = "#43aa8b", linewidth = 0.5) +  
  geom_line(data = plot_data_centile, aes(x = Age, y = P10), 
            color = "#90be6d", linewidth = 0.5) + 
  geom_line(data = plot_data_centile, aes(x = Age, y = P50), 
            color = "#f9c74f", linewidth = 0.5) +  
  geom_line(data = plot_data_centile, aes(x = Age, y = P90), 
            color = "#f8961e", linewidth = 0.5) +  
  geom_line(data = plot_data_centile, aes(x = Age, y = P95), 
            color = "#f3722c", linewidth = 0.5) +  
  labs(title = "Age-Normative Percentile Curves for Serum NfL",
       x = "Age (years)",
       y = "Serum NfL Level (pg/mL)",
       color = "Sex") +                                   # Change legend title
  scale_shape_manual(name = "Sex",
                     labels = c("Male", "Female"),
                     values = c(17, 16)) +  # 16 = Circle, 17 = Triangle
  scale_color_manual(                                     # Change legend categories
    values = c("0" = "#3D65A5", "1" = "#F05039"),              # Set custom colors
    labels = c("Male", "Female")                         # Change category labels
  ) +
  # Add labels for each curve
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)),  # Add labels at the max Age
            aes(x = Age, y = P97.5, label = "97.5"), 
            hjust = -0.2, color = "#f94144",size=2.5) +
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)), 
            aes(x = Age, y = P2.5, label = "2.5"), 
            hjust = -0.2, color = "#577590",size=2.5) +
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)), 
            aes(x = Age, y = P5, label = "5"), 
            hjust = -0.3, color = "#43aa8b",size=2.5) +
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)), 
            aes(x = Age, y = P10, label = "10"), 
            hjust = -0.2, color = "#90be6d",size=2.5) +
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)), 
            aes(x = Age, y = P50, label = "50"), 
            hjust = -0.2, color = "#f9c74f",size=2.5) +
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)), 
            aes(x = Age, y = P90, label = "90"), 
            hjust = -0.35, color = "#f8961e",size=2.5) +
  geom_text(data = plot_data_centile %>% filter(Age == max(Age)), 
            aes(x = Age, y = P95, label = "95"), 
            hjust = -0.2, color = "#f3722c",size=2.5) +
  theme_minimal()

#save --> graphical abstract fig 2
ggexport(P_by_age2, width=3800 ,height=3000, filename = "Graphic_ab_percentile_curves2.tiff", res = 550)

###---END---###