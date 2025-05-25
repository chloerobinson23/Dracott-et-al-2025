
## Foraging vs non foraging F-POD GAM models for vessel types
# 18 Oct 2024
# Author: Dr. Chloe V. Robinson (Based on analysis from Nuutilla et al. 2017)

RStudio.Version()
R.Version()

library(mgcv) #bam() function for models
library(car) #vif() function
library(itsadug)
library(tidymv)
library(dplyr)
library(ggplot2)

#read in data

hp <- read.csv(("hp_model_data.csv"), stringsAsFactors =TRUE, na.strings = c(""))
attach(hp)
str(hp)

#########
### Feeding minutes
########

#VIF testing

testmodel1<- lm(Feeding_Minutes~Month_Numeric + Hour + Vessel_Tanker + Vessel_Cargo + Vessel_Fishing + Vessel_Law
                + Vessel_Military + Vessel_Other + Vessel_Passenger + Vessel_Pilot + Vessel_Pleasure + Vessel_SAR
                + Vessel_Towing + Vessel_Towing.Large + Vessel_Tug + Vessel_Unknown + Vessel_Anti.pollution +
                  Speed_Slow + Speed_Med, data = hp)
summary(testmodel1)

#test vif
vif(testmodel1)

                       #GVIF^(1/(2*Df))
#Month_Numeric         1.010450  
#Hour                  1.067406
#Vessel_Tanker         1.000277 
#Vessel_Cargo          1.046757
#Vessel_Fishing        1.036317 
#Vessel_Law            1.080140
#Vessel_Military       1.000614
#Vessel_Other          1.008421
#Vessel_Passenger      1.015909
#Vessel_Pilot          1.100885
#Vessel_Pleasure       1.007539 
#Vessel_SAR            1.018972
#Vessel_Towing         1.019853
#Vessel_Towing.Large   1.001715
#Vessel_Tug            1.126719 
#Vessel_Unknown        1.003062  
#Vessel_Anti.pollution 1.011414
#Speed_Slow            1.014458 
#Speed_Med             1.247443 
#Speed_Fast            1.169629

## Use GVIF^(1/2) to make comparable to other vif estimates
## none above 5 so nothing to drop

### TESTS FOR AUTOCORRELATION

#run a base model - without autocorrelation

# removed vessel_tanker, vessel_towing.large, vessel_SAR, vessel_anti.pollution, vessel_military &
# vessel_unknown due to too low numbers to run for k

hpbambase2<- bam(Feeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24)
                +s(Vessel_Cargo, bs="tp", k=4) +s(Vessel_Fishing, bs="tp", k=3) +s(Vessel_Tug, bs="tp", k=5)
                +s(Vessel_Towing, bs="tp", k=2) +s(Vessel_Pilot, bs="tp", k=3)
                +s(Vessel_Pleasure, bs="tp", k=4) +s(Vessel_Passenger, bs="tp", k=3) 
                +s(Vessel_Other, bs="tp", k=3) +s(Speed_Slow, bs="tp", k=4) +s(Speed_Med, bs="tp", k=7)
                +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1))
#inspect correlation
par(mfrow=c(1,1), cex=1.1)


# default ACF function:
acf(resid(hpbambase2), main="acf(resid(hpbambase2))")

#Determine the value of lag 1, as indicated by the red dot in the picture below:
# we also ask to plot the ACF by specifying plot (FALSE by default):
r1hp <- start_value_rho(hpbambase2, plot=TRUE)

#this gives the r1 value
acf(resid(hpbambase2), plot=FALSE)$acf[2]
#r1 = 0.3702241

## above 0.2 so need to adjust


#determine Ar.start value
hp <- start_event(hp, column="Hour", event=c("Year","Date"),label.event="Event")

HP_new_Ar2 <-bam(Feeding_Minutes~s(Month_Numeric, bs="cs", k=12) + s(Hour, bs="cs",k=24)
                +s(Vessel_Cargo, bs="cs", k=4) +s(Vessel_Fishing, bs="cs", k=3) +s(Vessel_Tug, bs="cs", k=5)
                +s(Vessel_Towing, bs="cs", k=2) +s(Vessel_Pilot, bs="cs", k=3)
                +s(Vessel_Pleasure, bs="cs", k=4) +s(Vessel_Passenger, bs="cs", k=3) 
                +s(Vessel_Other, bs="cs", k=3) +s(Speed_Slow, bs="cs", k=4) +s(Speed_Med, bs="cs", k=7)
                +s(Speed_Fast, bs="cs", k=5) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

#To retrieve the corrected residuals of the model, one could use the function resid_gam.
par(mfrow=c(1,2), cex=1.1)

# normal residuals:
normal.res <- resid(HP_new_Ar2)
acf(normal.res, main="HP normal")

# corrected residuals:
corrected.res <- resid_gam(HP_new_Ar2)
acf(corrected.res,main="HP corrected")

## plot
summary(HP_new_Ar2)
plot(HP_new_Ar2, pages=1,all.terms=TRUE,ylim=c(-2, 2), shade=T, cex.axis=0.5, cex.lab=1)

#dealing with autocorrelation using bam - and negbin family instead of nb and theta value as calculated 
#by gam with nb - see the value in brackets in the gam summary (so run first a gam without AR structure,
#just to get the theta value)



##### MODELLING TIME!

#model1

model1<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24)
            +s(Vessel_Cargo, bs="tp", k=4) +s(Vessel_Fishing, bs="tp", k=3) +s(Vessel_Tug, bs="tp", k=5)
            +s(Vessel_Towing, bs="tp", k=2) +s(Vessel_Pilot, bs="tp", k=3)
            +s(Vessel_Pleasure, bs="tp", k=4) +s(Vessel_Passenger, bs="tp", k=3) 
            +s(Vessel_Other, bs="tp", k=3) +s(Speed_Slow, bs="tp", k=4) +s(Speed_Med, bs="tp", k=7)
            +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)


summary(model1)

#Parametric coefficients:
                                            #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         -0.65729    0.29502  -2.228   0.0259 *  
#as.factor(Year)2021 -0.37932    0.07282  -5.209 1.93e-07 ***
#as.factor(Year)2022 -0.98250    0.11655  -8.430  < 2e-16 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                  #edf Ref.df      F  p-value    
#s(Month_Numeric)    8.131 10.000 25.948  < 2e-16 ***
#s(Hour)             9.651 22.000 25.105  < 2e-16 ***
#s(Vessel_Cargo)     1.000  1.000  0.764 0.381973    
#s(Vessel_Fishing)   1.000  1.000  0.011 0.915537    
#s(Vessel_Tug)       2.601  3.017  6.504 0.000203 ***
#s(Vessel_Towing)    1.000  1.000  0.490 0.483903    
#s(Vessel_Pilot)     1.649  1.877  0.801 0.353224    
#s(Vessel_Pleasure)  1.000  1.000  0.075 0.784068    
#s(Vessel_Passenger) 1.000  1.000  2.446 0.117982    
#s(Vessel_Other)     1.000  1.000  4.683 0.030531 *  
#s(Speed_Slow)       1.000  1.000  0.555 0.456381    
#s(Speed_Med)        3.750  4.321 10.338  < 2e-16 ***
#s(Speed_Fast)       1.000  1.000 33.409  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.201   Deviance explained = 22.4%
#fREML =  11913  Scale est. = 0.63482   n = 12754


#model2 (remove non-sig values)

model2<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24) +s(Vessel_Tug, bs="tp", k=5)
            +s(Vessel_Other, bs="tp", k=3) +s(Speed_Med, bs="tp", k=7)  +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model2)

#Parametric coefficients:
                        #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          -1.02207    0.36316  -2.814  0.00489 **   
#as.factor(Year)2021 -0.38042    0.07256  -5.243 1.61e-07 ***
#as.factor(Year)2022 -0.98524    0.11625  -8.475  < 2e-16 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                #edf Ref.df      F  p-value    
#s(Month_Numeric) 8.121 10.000 25.850  < 2e-16 ***
#s(Hour)          9.662 22.000 25.399  < 2e-16 ***
#s(Vessel_Tug)    2.461  2.880  6.254 0.000321 ***
#s(Vessel_Other)  1.000  1.000  4.418 0.035623 *  
#s(Speed_Med)     3.870  4.432 14.888  < 2e-16 ***
#s(Speed_Fast)    1.130  1.246 15.372 1.75e-05 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =    0.2   Deviance explained = 22.3%
#fREML =  11898  Scale est. = 0.6321    n = 12754


#run AICs between 2 models
AIC(model1,model2)

#df      AIC
#model1 40.75033 36777.85
#model2 33.75141 36769.56


#plotting final HP DPM model
plot(model2,pages=1,all.terms=TRUE,ylim=c(-2, 2), shade=T,
     cex.axis=0.9, cex.lab=1.5, se=TRUE, shift = coef(model2)[1])


##### Plot the results of a Generalized Additive Model (GAM) on the response scale 
##### (rather than on the scale of the linear predictor or partial effects)####

### Month

# Check the structure of the data to verify factors are correctly set up
str(hp)

# Ensure all variables used in 's()' are numeric
hp$Month_Numeric <- as.numeric(hp$Month_Numeric)
hp$Hour <- as.numeric(hp$Hour)
hp$Vessel_Tug <- as.numeric(hp$Vessel_Tug)  # Ensure Vessel_Tug is numeric
hp$Vessel_Other <- as.numeric(hp$Vessel_Other)  # Ensure Vessel_Tug is numeric
hp$Speed_Med <- as.numeric(hp$Speed_Med)    # Ensure Speed_Med is numeric
hp$Speed_Fast <- as.numeric(hp$Speed_Fast)  # Ensure Speed_Fast is numeric

# Ensure 'Year' is a factor (for non-smoothed categorical effects)
hp$Year <- as.factor(hp$Year)

# Ensure 'Year' factor levels in new_data match those in the original model data
new_data$Year <- factor(new_data$Year, levels = levels(hp$Year))

# Drop any unused levels in 'Year' or any other factors
hp <- droplevels(hp)

####### ALL VARIABLES TOGETHER IN FACET ########

# Refit the model after ensuring the factors are correctly defined
model3 <- bam(Feeding_Minutes ~ 
                s(Month_Numeric, bs = "cc", k = 12) + 
                s(Hour, bs = "cc", k = 24) + 
                s(Vessel_Tug, bs = "tp", k = 5) +
                s(Vessel_Other, bs = "tp", k = 3) +
                s(Speed_Med, bs = "tp", k = 5) + 
                s(Speed_Fast, bs = "tp", k = 5), 
              data = hp, 
              family = negbin(0.1), 
              discrete = TRUE, 
              rho = r1hp, 
              AR.start = hp$start.event)

summary(model3)

# Create a data frame to store predictions for all variables

# Generate predictions for Month
new_data_month <- data.frame(
  Month_Numeric = seq(min(hp$Month_Numeric), max(hp$Month_Numeric), length.out = 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Tug = rep(mean(hp$Vessel_Tug), 300),
  Vessel_Other = rep(mean(hp$Vessel_Other), 300),
  Speed_Med = rep(mean(hp$Speed_Med), 300),
  Speed_Fast = rep(mean(hp$Speed_Fast), 300))

pred_month <- predict(model3, newdata = new_data_month, type = "response", se.fit = TRUE)
fit_month <- pred_month$fit
se_month <- pred_month$se.fit
lower_ci_month <- fit_month - 1.96 * se_month
upper_ci_month <- fit_month + 1.96 * se_month

# Store predictions for Month
df_month <- data.frame(
  Variable = "Month",
  Value = new_data_month$Month_Numeric,
  Fit = fit_month,
  Lower_CI = lower_ci_month,
  Upper_CI = upper_ci_month
)

# Generate predictions for Hour
new_data_hour <- data.frame(
  Month_Numeric = rep(mean(hp$Month_Numeric), 300),
  Hour = seq(min(hp$Hour), max(hp$Hour), length.out = 300),
  Vessel_Tug = rep(mean(hp$Vessel_Tug), 300),
  Vessel_Other = rep(mean(hp$Vessel_Other), 300),
  Speed_Med = rep(mean(hp$Speed_Med), 300),
  Speed_Fast = rep(mean(hp$Speed_Fast), 300)
)

pred_hour <- predict(model3, newdata = new_data_hour, type = "response", se.fit = TRUE)
fit_hour <- pred_hour$fit
se_hour <- pred_hour$se.fit
lower_ci_hour <- fit_hour - 1.96 * se_hour
upper_ci_hour <- fit_hour + 1.96 * se_hour

# Store predictions for Hour
df_hour <- data.frame(
  Variable = "Hour",
  Value = new_data_hour$Hour,
  Fit = fit_hour,
  Lower_CI = lower_ci_hour,
  Upper_CI = upper_ci_hour
)

# Generate predictions for Vessel_Tug
new_data_tug <- data.frame(
  Month_Numeric = rep(mean(hp$Month_Numeric), 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Tug = seq(min(hp$Vessel_Tug), max(hp$Vessel_Tug), length.out = 300),
  Vessel_Other = rep(mean(hp$Vessel_Other), 300),
  Speed_Med = rep(mean(hp$Speed_Med), 300),
  Speed_Fast = rep(mean(hp$Speed_Fast), 300)
)

pred_tug <- predict(model3, newdata = new_data_tug, type = "response", se.fit = TRUE)
fit_tug <- pred_tug$fit
se_tug <- pred_tug$se.fit
lower_ci_tug <- fit_tug - 1.96 * se_tug
upper_ci_tug <- fit_tug + 1.96 * se_tug

# Store predictions for Vessel_Tug
df_tug <- data.frame(
  Variable = "Vessel_Tug",
  Value = new_data_tug$Vessel_Tug,
  Fit = fit_tug,
  Lower_CI = lower_ci_tug,
  Upper_CI = upper_ci_tug
)

# Generate predictions for Vessel_Other
new_data_other <- data.frame(
  Month_Numeric = rep(mean(hp$Month_Numeric), 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Tug = rep(mean(hp$Vessel_Tug), 300),
  Vessel_Other = seq(min(hp$Speed_Med), max(hp$Speed_Med),length.out =300),
  Speed_Med = rep(mean(hp$Speed_Med), 300),
  Speed_Fast = rep(mean(hp$Speed_Fast), 300)
)

pred_other <- predict(model3, newdata = new_data_other, type = "response", se.fit = TRUE)
fit_other <- pred_other$fit
se_other <- pred_other$se.fit
lower_ci_other <- fit_other - 1.96 * se_tug
upper_ci_other <- fit_other + 1.96 * se_tug

# Store predictions for Vessel_Other
df_other <- data.frame(
  Variable = "Vessel_Other",
  Value = new_data_other$Vessel_Other,
  Fit = fit_other,
  Lower_CI = lower_ci_other,
  Upper_CI = upper_ci_other
)

# Generate predictions for Speed_Med
new_data_speed_med <- data.frame(
  Month_Numeric = rep(mean(hp$Month_Numeric), 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Tug = rep(mean(hp$Vessel_Tug), 300),
  Vessel_Other = rep(mean(hp$Vessel_Other), 300),
  Speed_Med = seq(min(hp$Speed_Med), max(hp$Speed_Med), length.out = 300),
  Speed_Fast = rep(mean(hp$Speed_Fast), 300)
)

pred_speed_med <- predict(model3, newdata = new_data_speed_med, type = "response", se.fit = TRUE)
fit_speed_med <- pred_speed_med$fit
se_speed_med <- pred_speed_med$se.fit
lower_ci_speed_med <- fit_speed_med - 1.96 * se_speed_med
upper_ci_speed_med <- fit_speed_med + 1.96 * se_speed_med

# Store predictions for Speed_Med
df_speed_med <- data.frame(
  Variable = "Speed_Med",
  Value = new_data_speed_med$Speed_Med,
  Fit = fit_speed_med,
  Lower_CI = lower_ci_speed_med,
  Upper_CI = upper_ci_speed_med
)

# Generate predictions for Speed_Fast
new_data_speed_fast <- data.frame(
  Month_Numeric = rep(mean(hp$Month_Numeric), 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Tug = rep(mean(hp$Vessel_Tug), 300),
  Vessel_Other = rep(mean(hp$Vessel_Other), 300),
  Speed_Med = rep(mean(hp$Speed_Med), 300),
  Speed_Fast = seq(min(hp$Speed_Fast), max(hp$Speed_Fast), length.out = 300)
)

pred_speed_fast <- predict(model3, newdata = new_data_speed_fast, type = "response", se.fit = TRUE)
fit_speed_fast <- pred_speed_fast$fit
se_speed_fast <- pred_speed_fast$se.fit
lower_ci_speed_fast <- fit_speed_fast - 1.96 * se_speed_fast
upper_ci_speed_fast <- fit_speed_fast + 1.96 * se_speed_fast

# Store predictions for Speed_Fast
df_speed_fast <- data.frame(
  Variable = "Speed_Fast",
  Value = new_data_speed_fast$Speed_Fast,
  Fit = fit_speed_fast,
  Lower_CI = lower_ci_speed_fast,
  Upper_CI = upper_ci_speed_fast
)

# Combine all data frames into one for ggplot
df_all <- rbind(df_month, df_hour, df_tug, df_other, df_speed_med, df_speed_fast)


# Reorder 'Variable' in the dataframe
df_all$Variable <- factor(df_all$Variable, levels = c("Month", "Hour", "Vessel_Tug", "Vessel_Other", "Speed_Med", "Speed_Fast"))

# Plot using ggplot2 with faceting
figure1 <- ggplot(df_all, aes(x = Value, y = Fit)) +
  geom_line(color = "blue", size = 1) +  # Use geom_line instead of geom_smooth
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), 
              fill = "deepskyblue3", alpha = 0.3) + 
  facet_wrap(~ Variable, scales = "free_x") +  # Facet by reordered variable
  labs(x = "Value", y = "Predicted Foraging Minutes", 
       title = "") +
  theme_minimal() +  # Minimal theme
  theme(strip.text = element_text(size = 12),  # Smaller facet labels
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 14, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14),  # Center and size for title
        legend.position = "none",  # Remove legend
        panel.grid = element_blank(),  # Remove gridlines
        panel.background = element_blank(),  # Remove the panel background
        strip.background = element_blank(),  # Remove facet background
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add borders around each facet
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

figure1

#save plot
ggsave("Figure1_Foraging_DPM_Vessel.tiff", figure1, dpi=400)



##difference in explained deviance
1-(model1$dev/model1$null)
#0.2241069

##difference in explained deviance
1-(model2$dev/model2$null)
#0.2232133


#perform anova to compare models (start with model with least variables)
anova(model2, model1, test="Chisq")

#Resid. Df Resid. Dev       Df Deviance  Pr(>Chi)    
#1     12716     4960.4                         
#2     12709     4954.6 7.1115   5.7059   0.2631


### accept null hypothesis that all models are the same 

########################################


########################
## AIC table 
#######################

#run model 4 as if carrying on from model 3, removing the vars that we know are signif
#remove month first

model4<-bam(Feeding_Minutes~s(Hour, bs="cc",k=24) +s(Vessel_Tug, bs="tp", k=5)
            +s(Vessel_Other, bs="tp", k=3) +s(Speed_Med, bs="tp", k=7) +s(Speed_Fast, bs="tp", k=5)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model4)

#now remove hour and add month back in 

model5<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc",k=12) +s(Vessel_Tug, bs="tp", k=5)
            +s(Vessel_Other, bs="tp", k=3) +s(Speed_Med, bs="tp", k=7) +s(Speed_Fast, bs="tp", k=5)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model5)

#now remove vessel_tug and add hour back in 

model6<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc",k=12) +s(Hour, bs="cc",k=24)
            +s(Vessel_Other, bs="tp", k=3) +s(Speed_Med, bs="tp", k=7) +s(Speed_Fast, bs="tp", k=5)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model6)

#now remove vessel_other and add vessel_tug back in

model7<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc",k=12) +s(Hour, bs="cc",k=24)
            +s(Vessel_Tug, bs="tp", k=5) +s(Speed_Med, bs="tp", k=7) +s(Speed_Fast, bs="tp", k=5)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model7)

#now remove speed_med and add vessel other back in 

model8<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc",k=12) +s(Hour, bs="cc",k=24)
            +s(Vessel_Tug, bs="tp", k=5) +s(Vessel_Other, bs="tp", k=3) +s(Speed_Fast, bs="tp", k=5)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model8)

#now remove speed_fast and add speed_med back in 

model9<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc",k=12) +s(Hour, bs="cc",k=24)
            +s(Vessel_Tug, bs="tp", k=5) +s(Vessel_Other, bs="tp", k=3) +s(Speed_Med, bs="tp", k=7)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model9)

#now remove year and add speed_fast back in 

model10<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc",k=12) +s(Hour, bs="cc",k=24)
            +s(Vessel_Tug, bs="tp", k=5) +s(Vessel_Other, bs="tp", k=3) +s(Speed_Med, bs="tp", k=7)
            +s(Speed_Fast, bs="tp", k=5), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model10)



#run AICs
AIC(model2, model4,model5, model6,model7, model8, model9, model10)
#df   AIC 
#model2  33.75141 36769.56
#model4  27.15007 37036.38
#model5  23.21570 37330.67
#model6  30.44825 36775.22
#model7  32.67456 36771.38
#model8  29.16170 36823.73
#model9  31.04156 36804.11
#model10 34.53012 36850.75

##################################




###########################
### Non-feeding minutes
############################

#VIF testing

testmodel2<- lm(NonFeeding_Minutes~Month_Numeric + Hour + Vessel_Tanker + Vessel_Cargo + Vessel_Fishing + Vessel_Law
                + Vessel_Military + Vessel_Other + Vessel_Passenger + Vessel_Pilot + Vessel_Pleasure + Vessel_SAR
                + Vessel_Towing + Vessel_Towing.Large + Vessel_Tug + Vessel_Unknown + Vessel_Anti.pollution +
                  Speed_Slow + Speed_Med +Speed_Fast, data = hp)
summary(testmodel2)

#test vif
vif(testmodel2)


## Use GVIF^(1/2) to make comparable to other vif estimates
## none above 5 so nothing to drop

### TESTS FOR AUTOCORRELATION

#run a base model - without autocorrelation

# removed vessel_tanker, vessel_towing.large, vessel_SAR, vessel_anti.pollution, vessel_military &
# vessel_unknown due to too low numbers to run for k

hpbambase2<- bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24)
                 +s(Vessel_Cargo, bs="tp", k=4) +s(Vessel_Fishing, bs="tp", k=3) +s(Vessel_Tug, bs="tp", k=5)
                 +s(Vessel_Towing, bs="tp", k=2) +s(Vessel_Pilot, bs="tp", k=3)
                 +s(Vessel_Pleasure, bs="tp", k=4) +s(Vessel_Passenger, bs="tp", k=3) 
                 +s(Vessel_Other, bs="tp", k=3) +s(Speed_Slow, bs="tp", k=4) +s(Speed_Med, bs="tp", k=7)
                 +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1))
#inspect correlation
par(mfrow=c(1,1), cex=1.1)


# default ACF function:
acf(resid(hpbambase2), main="acf(resid(hpbambase2))")

#Determine the value of lag 1, as indicated by the red dot in the picture below:
# we also ask to plot the ACF by specifying plot (FALSE by default):
r1hp <- start_value_rho(hpbambase2, plot=TRUE)

#this gives the r1 value
acf(resid(hpbambase2), plot=FALSE)$acf[2]
#r1 = 0.4142446

## above 0.2 so need to adjust


#determine Ar.start value
hp <- start_event(hp, column="Hour", event=c("Year","Date"),label.event="Event")

HP_new_Ar2 <-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24)
                 +s(Vessel_Cargo, bs="tp", k=4) +s(Vessel_Fishing, bs="tp", k=3) +s(Vessel_Tug, bs="tp", k=5)
                 +s(Vessel_Towing, bs="tp", k=2) +s(Vessel_Pilot, bs="tp", k=3)
                 +s(Vessel_Pleasure, bs="tp", k=4) +s(Vessel_Passenger, bs="tp", k=3) 
                 +s(Vessel_Other, bs="tp", k=3) +s(Speed_Slow, bs="tp", k=4) +s(Speed_Med, bs="tp", k=7)
                 +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

#To retrieve the corrected residuals of the model, one could use the function resid_gam.
par(mfrow=c(1,2), cex=1.1)

# normal residuals:
normal.res <- resid(HP_new_Ar2)
acf(normal.res, main="HP normal")

# corrected residuals:
corrected.res <- resid_gam(HP_new_Ar2)
acf(corrected.res,main="HP corrected")

#dealing with autocorrelation using bam - and negbin family instead of nb and theta value as calculated 
#by gam with nb - see the value in brackets in the gam summary (so run first a gam without AR structure,
#just to get the theta value)



##### MODELLING TIME!

#model4

model4<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24)
            +s(Vessel_Cargo, bs="tp", k=4) +s(Vessel_Fishing, bs="tp", k=3) +s(Vessel_Tug, bs="tp", k=5)
            +s(Vessel_Towing, bs="tp", k=3) +s(Vessel_Pilot, bs="tp", k=3)
            +s(Vessel_Pleasure, bs="tp", k=4) +s(Vessel_Passenger, bs="tp", k=3) 
            +s(Vessel_Other, bs="tp", k=3) +s(Speed_Slow, bs="tp", k=4) +s(Speed_Med, bs="tp", k=7)
            +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)


summary(model4)

#Parametric coefficients:
                    #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         -0.67887    0.42407  -1.601    0.109    
#as.factor(Year)2021 -0.25214    0.05248  -4.804 1.57e-06 ***
#as.factor(Year)2022 -0.51428    0.08147  -6.312 2.84e-10 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                    #edf Ref.df      F  p-value    
#s(Month_Numeric)    8.419 10.000 33.887  < 2e-16 ***
#s(Hour)             9.331 22.000 17.474  < 2e-16 ***
#s(Vessel_Cargo)     1.000  1.001  4.208  0.04000 *  
#s(Vessel_Fishing)   1.000  1.000  1.237  0.26657    
#s(Vessel_Tug)       2.189  2.604  9.317 5.45e-05 ***
#s(Vessel_Towing)    1.010  1.020  1.544  0.22781    
#s(Vessel_Pilot)     1.759  1.942  3.103  0.03351 *  
#s(Vessel_Pleasure)  1.000  1.000  3.546  0.05971 .  
#s(Vessel_Passenger) 1.000  1.000 12.101  0.00051 ***
#s(Vessel_Other)     1.000  1.000  1.224  0.26892    
#s(Speed_Slow)       1.000  1.000  1.729  0.18876    
#s(Speed_Med)        3.882  4.443 19.160  < 2e-16 ***
#s(Speed_Fast)       1.928  2.303 31.482  < 2e-16 ***
  #---

#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.195   Deviance explained = 15.8%
#fREML = 9621.8  Scale est. = 0.32076   n = 12754


#model5 (remove non-sig values)

model5<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24)
            +s(Vessel_Cargo, bs="tp", k=4)  +s(Vessel_Tug, bs="tp", k=5) +s(Vessel_Pilot, bs="tp", k=3) +s(Vessel_Passenger, bs="tp", k=3) +s(Speed_Med, bs="tp", k=7)
            +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

summary(model5)

#Parametric coefficients:
                    #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.22781    0.24730   0.921    0.357    
#as.factor(Year)2021 -0.25714    0.05196  -4.949 7.57e-07 ***
#as.factor(Year)2022 -0.51736    0.08075  -6.407 1.53e-10 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                #edf Ref.df      F  p-value    
#s(Month_Numeric)    8.436 10.000 34.792  < 2e-16 ***
#s(Hour)             9.362 22.000 17.844  < 2e-16 ***
#s(Vessel_Cargo)     1.000  1.000  3.907   0.0480 *  
#s(Vessel_Tug)       2.112  2.522  9.267 6.23e-05 ***
#s(Vessel_Pilot)     1.715  1.918  2.404   0.0632 .  
#s(Vessel_Passenger) 1.000  1.000 10.661   0.0011 ** 
#s(Speed_Med)        3.947  4.504 23.791  < 2e-16 ***
#s(Speed_Fast)       2.019  2.401 38.222  < 2e-16 ***
  #---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.195   Deviance explained = 15.8%
#fREML = 9654.1  Scale est. = 0.31897   n = 12754


#model6 (remove non-sig values)

model6<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24)
            +s(Vessel_Tug, bs="tp", k=5) +s(Vessel_Passenger, bs="tp", k=3) 
            +s(Speed_Med, bs="tp", k=7) +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

summary(model6)

#Parametric coefficients:
                      #Estimate Std. Error t value Pr(>|t|)
#(Intercept)          0.35580    0.23750   1.498    0.134    
#as.factor(Year)2021 -0.25578    0.05247  -4.875 1.10e-06 ***
#as.factor(Year)2022 -0.51655    0.08152  -6.337 2.43e-10 ***
  #---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                #edf Ref.df      F  p-value    
#s(Month_Numeric)    8.406 10.000 34.186  < 2e-16 ***
#s(Hour)             9.410 22.000 17.461  < 2e-16 ***
#s(Vessel_Tug)       1.976  2.377  8.348 0.000263 ***
#s(Vessel_Passenger) 1.000  1.000  8.790 0.003057 ** 
#s(Speed_Med)        4.048  4.597 31.835  < 2e-16 ***
#s(Speed_Fast)       2.171  2.566 46.464  < 2e-16 ***
  #---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.195   Deviance explained = 15.7%
#fREML = 9621.8  Scale est. = 0.32137   n = 12754


#run AICs between 2 models
AIC(model4,model5,model6)

#df      AIC
#model4 40.97050 64703.89
#model5 36.05824 64767.08
#model6 32.94775 64695.17


#AIC for model6 is lower (and contains only significant variables)

##difference in explained deviance
1-(model6$dev/model6$null)
#0.1566114

##difference in explained deviance
1-(model5$dev/model5$null)
#0.157649

##difference in explained deviance
1-(model4$dev/model4$null)
#0.1581615

#perform anova to compare models (start with model with least variables)
anova(model6, model5, model4, test="Chisq")

#Resid. Df Resid. Dev       Df Deviance  Pr(>Chi)    
#1     12717     3984.8                            
#2     12713     3979.9 3.2908   4.9022 0.002144 **
#3     12709     3977.4 4.9099   2.4216 0.175498   


### reject null hypothesis that all models are the same 


####### ALL VARIABLES TOGETHER IN FACET ########

# Refit the model after ensuring the factors are correctly defined (remove year for plotting smooth only)
model6<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24)
            +s(Vessel_Tug, bs="tp", k=5) +s(Vessel_Passenger, bs="tp", k=3) 
            +s(Speed_Med, bs="tp", k=7) +s(Speed_Fast, bs="tp", k=5), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

# Create a data frame to store predictions for all variables

# Generate predictions for Month
new_data_month <- data.frame(
  Month_Numeric = seq(min(hp$Month_Numeric), max(hp$Month_Numeric), length.out = 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Cargo = rep(mean(hp$Vessel_Other), 300),
  Vessel_Tug = rep(mean(hp$Vessel_Tug), 300),
  Vessel_Passenger = rep(mean(hp$Vessel_Other), 300),
  Speed_Med = rep(mean(hp$Speed_Med), 300),
  Speed_Fast = rep(mean(hp$Speed_Fast), 300))

pred_month <- predict(model6, newdata = new_data_month, type = "response", se.fit = TRUE)
fit_month <- pred_month$fit
se_month <- pred_month$se.fit
lower_ci_month <- fit_month - 1.96 * se_month
upper_ci_month <- fit_month + 1.96 * se_month

# Store predictions for Month
df_month <- data.frame(
  Variable = "Month",
  Value = new_data_month$Month_Numeric,
  Fit = fit_month,
  Lower_CI = lower_ci_month,
  Upper_CI = upper_ci_month
)

# Generate predictions for Hour
new_data_hour <- data.frame(
  Month_Numeric = rep(mean(hp$Month_Numeric), 300),
  Hour = seq(min(hp$Hour), max(hp$Hour), length.out = 300),
  Vessel_Cargo = rep(mean(hp$Vessel_Other), 300),
  Vessel_Tug = rep(mean(hp$Vessel_Tug), 300),
  Vessel_Passenger = rep(mean(hp$Vessel_Other), 300),
  Speed_Med = rep(mean(hp$Speed_Med), 300),
  Speed_Fast = rep(mean(hp$Speed_Fast), 300))

pred_hour <- predict(model6, newdata = new_data_hour, type = "response", se.fit = TRUE)
fit_hour <- pred_hour$fit
se_hour <- pred_hour$se.fit
lower_ci_hour <- fit_hour - 1.96 * se_hour
upper_ci_hour <- fit_hour + 1.96 * se_hour

# Store predictions for Hour
df_hour <- data.frame(
  Variable = "Hour",
  Value = new_data_hour$Hour,
  Fit = fit_hour,
  Lower_CI = lower_ci_hour,
  Upper_CI = upper_ci_hour
)


# Generate predictions for Vessel_Tug
new_data_vessel_tug <- data.frame(
  Month_Numeric = rep(mean(hp$Month_Numeric), 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Cargo = rep(mean(hp$Vessel_Cargo), 300),
  Vessel_Tug = seq(min(hp$Vessel_Tug), max(hp$Vessel_Tug), length.out = 300),
  Vessel_Passenger = rep(mean(hp$Vessel_Other), 300),
  Speed_Med = rep(mean(hp$Speed_Med), 300),
  Speed_Fast = rep(mean(hp$Speed_Fast), 300))

pred_vessel_tug <- predict(model6, newdata = new_data_vessel_tug, type = "response", se.fit = TRUE)
fit_vessel_tug <- pred_vessel_tug$fit
se_vessel_tug <- pred_vessel_tug$se.fit
lower_ci_vessel_tug <- fit_vessel_tug - 1.96 * se_vessel_tug
upper_ci_vessel_tug <- fit_vessel_tug + 1.96 * se_vessel_tug

# Store predictions for Vessel_Tug
df_vessel_tug <- data.frame(
  Variable = "Vessel_Tug",
  Value = new_data_vessel_tug$Vessel_Tug,
  Fit = fit_vessel_tug,
  Lower_CI = lower_ci_vessel_tug,
  Upper_CI = upper_ci_vessel_tug
)

# Generate predictions for Vessel_Passenger
new_data_vessel_passenger <- data.frame(
  Month_Numeric = rep(mean(hp$Month_Numeric), 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Cargo = rep(mean(hp$Vessel_Cargo), 300),
  Vessel_Tug = rep(mean(hp$Vessel_Tug), 300),
  Vessel_Passenger = seq(min(hp$Vessel_Other), max(hp$Vessel_Passenger), length.out = 300),
  Speed_Med = rep(mean(hp$Speed_Med), 300),
  Speed_Fast = rep(mean(hp$Speed_Fast), 300))

pred_vessel_passenger <- predict(model6, newdata = new_data_vessel_passenger, type = "response", se.fit = TRUE)
fit_vessel_passenger <- pred_vessel_passenger$fit
se_vessel_passenger <- pred_vessel_passenger$se.fit
lower_ci_vessel_passenger <- fit_vessel_passenger - 1.96 * se_vessel_passenger
upper_ci_vessel_passenger <- fit_vessel_passenger + 1.96 * se_vessel_passenger

# Store predictions for Vessel_Passenger
df_vessel_passenger <- data.frame(
  Variable = "Vessel_Passenger",
  Value = new_data_vessel_passenger$Vessel_Passenger,
  Fit = fit_vessel_passenger,
  Lower_CI = lower_ci_vessel_passenger,
  Upper_CI = upper_ci_vessel_passenger
)

# Generate predictions for Speed_Med
new_data_speed_med <- data.frame(
  Month_Numeric = rep(mean(hp$Month_Numeric), 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Cargo = rep(mean(hp$Vessel_Other), 300),
  Vessel_Tug = rep(mean(hp$Vessel_Tug), 300),
  Vessel_Passenger = rep(mean(hp$Vessel_Other), 300),
  Speed_Med = seq(min(hp$Speed_Med), max(hp$Speed_Med), length.out = 300),
  Speed_Fast = rep(mean(hp$Speed_Fast), 300))

pred_speed_med <- predict(model6, newdata = new_data_speed_med, type = "response", se.fit = TRUE)
fit_speed_med <- pred_speed_med$fit
se_speed_med <- pred_speed_med$se.fit
lower_ci_speed_med <- fit_speed_med - 1.96 * se_speed_med
upper_ci_speed_med <- fit_speed_med + 1.96 * se_speed_med

# Store predictions for Speed_Med
df_speed_med <- data.frame(
  Variable = "Speed_Med",
  Value = new_data_speed_med$Speed_Med,
  Fit = fit_speed_med,
  Lower_CI = lower_ci_speed_med,
  Upper_CI = upper_ci_speed_med
)

# Generate predictions for Speed_Fast
new_data_speed_fast <- data.frame(
  Month_Numeric = rep(mean(hp$Month_Numeric), 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Cargo = rep(mean(hp$Vessel_Other), 300),
  Vessel_Tug = rep(mean(hp$Vessel_Tug), 300),
  Vessel_Passenger = rep(mean(hp$Vessel_Other), 300),
  Speed_Med = rep(mean(hp$Speed_Med), 300),
  Speed_Fast = seq(min(hp$Speed_Fast), max(hp$Speed_Fast), length.out = 300))

pred_speed_fast <- predict(model6, newdata = new_data_speed_fast, type = "response", se.fit = TRUE)
fit_speed_fast <- pred_speed_fast$fit
se_speed_fast <- pred_speed_fast$se.fit
lower_ci_speed_fast <- fit_speed_fast - 1.96 * se_speed_fast
upper_ci_speed_fast <- fit_speed_fast + 1.96 * se_speed_fast

# Store predictions for Speed_Fast
df_speed_fast <- data.frame(
  Variable = "Speed_Fast",
  Value = new_data_speed_fast$Speed_Fast,
  Fit = fit_speed_fast,
  Lower_CI = lower_ci_speed_fast,
  Upper_CI = upper_ci_speed_fast
)

# Combine all data frames into one for ggplot
df_all <- rbind(df_month, df_hour, df_vessel_tug, df_vessel_passenger, df_speed_med, df_speed_fast)


# Reorder 'Variable' in the dataframe
df_all$Variable <- factor(df_all$Variable, levels = c("Month", "Hour", "Vessel_Tug", "Vessel_Passenger", "Speed_Med", "Speed_Fast"))

dev.off()

# Plot using ggplot2 with faceting
figure2 <-ggplot(df_all, aes(x = Value, y = Fit)) +
  geom_line(color = "green", size = 1) +  # Use geom_line instead of geom_smooth
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), 
              fill = "#90EE90", alpha = 0.3) + 
  facet_wrap(~ Variable, scales = "free_x") +  # Facet by variable
  labs(x = "Value", y = "Predicted Non-Foraging Minutes", 
       title = "") +
  theme_minimal() +  # Minimal theme
  theme(strip.text = element_text(size = 12),  # Smaller facet labels
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 14, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14),  # Center and size for title
        legend.position = "none",  # Remove legend
        panel.grid = element_blank(),  # Remove gridlines
        panel.background = element_blank(),  # Remove the panel background
        strip.background = element_blank(),  # Remove facet background
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add borders around each facet
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

figure2

#save plot
ggsave("Figure2_Non-Foraging_DPM_Vessel.tiff", figure2, dpi=400)


########################
## AIC table
#######################

#run model 7 as if carrying on from model 6, removing the vars that we know are signif
#remove month first

model7<-bam(NonFeeding_Minutes~s(Hour, bs="cc",k=24)
            +s(Vessel_Tug, bs="tp", k=5) +s(Vessel_Passenger, bs="tp", k=3) 
            +s(Speed_Med, bs="tp", k=7) +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model7)

#now remove hour and add month back in 

model8<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc",k=12)
            +s(Vessel_Tug, bs="tp", k=5) +s(Vessel_Passenger, bs="tp", k=3) 
            +s(Speed_Med, bs="tp", k=7) +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

summary(model8)


#now remove vessel_tug and add vessel_cargo back in

model9<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12)+s(Hour, bs="cc",k=24) 
            +s(Vessel_Passenger, bs="tp", k=3) 
            +s(Speed_Med, bs="tp", k=7) +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

summary(model9)

#now remove vessel_passenger and add vessel tug back in 


model10<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12)+s(Hour, bs="cc",k=24) 
             +s(Vessel_Tug, bs="tp", k=5) 
             +s(Speed_Med, bs="tp", k=7) +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

summary(model10)


#now remove speed_med and add vessel passenger back in 

model11<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12)+s(Hour, bs="cc",k=24) 
             +s(Vessel_Tug, bs="tp", k=5) 
             +s(Vessel_Passenger, bs="tp", k=3) +s(Speed_Fast, bs="tp", k=5) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

summary(model11)


#now remove speed_fast and add speed_med back in 

model12<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12)+s(Hour, bs="cc",k=24) 
             +s(Vessel_Tug, bs="tp", k=5) 
             +s(Vessel_Passenger, bs="tp", k=3) +s(Speed_Med, bs="tp", k=7) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

summary(model12)

#now remove year and add speed_fast back in 

model13<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12)+s(Hour, bs="cc",k=24) 
            +s(Vessel_Tug, bs="tp", k=5) 
             +s(Vessel_Passenger, bs="tp", k=3) +s(Speed_Med, bs="tp", k=7) +s(Speed_Fast, bs="tp", k=5), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

summary(model13)



#run AICs
AIC(model6, model7, model8,model9, model10, model11, model12, model13)
#df   AIC 
#model6  32.94775 64695.17
#model7  27.43033 64902.13
#model8  23.54500 64906.53
#model9  30.50450 64695.81
#model10 32.04269 64695.09
#model11 28.79489 64751.08
#model12 30.47761 64735.36
#model13 30.91796 64721.21

##################################
### END
#################################
