
## Foraging vs non foraging F-POD GAM models for vessel numbers
# 18 Oct 2024
# Author: Dr. Chloe V. Robinson (Based on analysis from Nuutilla et al. 2017)

RStudio.Version()
R.Version()

library(mgcv) #bam() function for models
library(car) #vif() function
library(itsadug)

#read in data

hp <- read.csv(("hp_model_data.csv"), stringsAsFactors =TRUE, na.strings = c(""))
attach(hp)
str(hp)

#########
### Feeding minutes
########

#VIF testing

testmodel1<- lm(Feeding_Minutes~Month_Numeric + Hour + Vessel_Total + Speed_Fast + Speed_Med + Speed_Slow, data = hp)
summary(testmodel1)

#test vif
vif(testmodel1)

                #GVIF Df GVIF^(1/(2*Df))
#Month_Numeric        1.033935 11        1.006169 
#Hour         1.010923  1        1.008505 
#Vessel_Total 1.594274  1        1.590953
#Speed_Fast   1.286605  1        1.279216 
#Speed_Med    1.343427  1        1.340535 
#Speed_Slow   1.025765  1        1.019583 

## Use GVIF^(1/2) to make comparable to other vif estimates
## none above 5 so nothing to drop

### TESTS FOR AUTOCORRELATION

#run a base model - without autocorrelation

hpbambase<- bam(Feeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24) +s(Vessel_Total, bs="tp", k=6) +s(Speed_Fast, bs="tp", k=5) 
                +s(Speed_Med, bs="tp",k=6) + s(Speed_Slow, bs="tp", k=4) +as.factor(Year), data=hp, family=negbin(0.1))
#inspect correlation
par(mfrow=c(1,1), cex=1.1)


# default ACF function:
acf(resid(hpbambase), main="acf(resid(hpbambase))")

#Determine the value of lag 1, as indicated by the red dot in the picture below:
# we also ask to plot the ACF by specifying plot (FALSE by default):
r1hp <- start_value_rho(hpbambase, plot=TRUE)

#this gives the r1 value
acf(resid(hpbambase), plot=FALSE)$acf[2]
#r1 =0.3706485

### above 0.2 so need to run more tests


#determine Ar.start value
hp <- start_event(hp, column="Hour", event=c("Year","Date"),label.event="Event")

HP_new_Ar <-bam(Feeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24) +s(Vessel_Total, bs="tp", k=6) +s(Speed_Fast, bs="tp", k=5) 
                +s(Speed_Med, bs="tp",k=6) + s(Speed_Slow, bs="tp", k=4) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

#To retrieve the corrected residuals of the model, one could use the function resid_gam.
par(mfrow=c(1,2), cex=1.1)

# normal residuals:
normal.res <- resid(HP_new_Ar)
acf(normal.res, main="HP normal")

# corrected residuals:
corrected.res <- resid_gam(HP_new_Ar)
acf(corrected.res,main="HP corrected")


##### MODELLING TIME!

#model1

model1<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24) +s(Vessel_Total, bs="tp", k=6) +s(Speed_Fast, bs="tp", k=5) 
                +s(Speed_Med, bs="tp",k=6) + s(Speed_Slow, bs="tp", k=4) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model1)

#Parametric coefficients:
                                #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.21373    0.27321   0.782    0.434    
#as.factor(Year)2021 -0.36238    0.07316  -4.953 7.39e-07 ***
#as.factor(Year)2022 -0.97816    0.11739  -8.332  < 2e-16 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                #edf Ref.df      F p-value    
#s(Month_Numeric) 8.097 10.000 25.395  <2e-16 ***
#s(Hour)          9.613 22.000 25.106  <2e-16 ***
#s(Vessel_Total)  3.906  4.389 17.311  <2e-16 ***
#s(Speed_Fast)    1.000  1.000  0.334   0.564    
#s(Speed_Med)     1.000  1.000  0.001   0.977    
#s(Speed_Slow)    1.000  1.001  1.346   0.246    
#---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.201   Deviance explained = 22.6%
#fREML =  11921  Scale est. = 0.64271   n = 12754


#model2 (remove speed_fast, speed_med and speed_slow)

model2<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24) +s(Vessel_Total, bs="tp", k=6)
         +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model2)

#Parametric coefficients:
                      #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.06979    0.09414   0.741    0.459    
#as.factor(Year)2021 -0.36332    0.07332  -4.955 7.32e-07 ***
#as.factor(Year)2022 -0.98002    0.11769  -8.327  < 2e-16 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                #edf Ref.df     F p-value    
#s(Month_Numeric) 8.085 10.000 24.97  <2e-16 ***
#s(Hour)          9.612 22.000 25.14  <2e-16 ***
#s(Vessel_Total)  3.912  4.385 41.64  <2e-16 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.201   Deviance explained = 22.6%
#fREML =  11929  Scale est. = 0.64594   n = 12754


#run AICs between 2 models
AIC(model1,model2)

            #df      AIC
#model1 29.42636 36739.11
#model2 26.38844 36732.19


#AIC for model2 marginally smaller and contains non-sig variables = better

##difference in explained deviance
1-(model1$dev/model1$null)
#0.2259477

##difference in explained deviance
1-(model2$dev/model2$null)
#0.2260808


#perform anova to compare models (start with model with least variables)
anova(model2, model1, test="Chisq")

      #Resid. Df Resid. Dev       Df Deviance  Pr(>Chi)    
#1     12724     4942.0                         
#2     12721     4942.9 3.0458 -0.84978  


### accept null hypothesis that all models are the same 

### export the plot

####### ALL VARIABLES TOGETHER IN FACET ########

# Refit the model after ensuring the factors are correctly defined
model3<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24) +s(Vessel_Total, bs="tp", k=6),
             data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

summary(model3)

# Create a data frame to store predictions for all variables

# Generate predictions for Month
new_data_month <- data.frame(
  Month_Numeric = seq(min(hp$Month_Numeric), max(hp$Month_Numeric), length.out = 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Total = rep(mean(hp$Vessel_Total), 300))

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
  Vessel_Total = rep(mean(hp$Vessel_Total), 300)
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

# Generate predictions for Vessel_Total
new_data_total <- data.frame(
  Month_Numeric = rep(mean(hp$Month_Numeric), 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Total = seq(min(hp$Vessel_Total), max(hp$Vessel_Total), length.out = 300)
)

pred_total <- predict(model3, newdata = new_data_total, type = "response", se.fit = TRUE)
fit_total <- pred_total$fit
se_total <- pred_total$se.fit
lower_ci_total <- fit_total - 1.96 * se_total
upper_ci_total <- fit_total + 1.96 * se_total

# Store predictions for Vessel_total
df_total <- data.frame(
  Variable = "Vessel_total",
  Value = new_data_total$Vessel_Total,
  Fit = fit_total,
  Lower_CI = lower_ci_total,
  Upper_CI = upper_ci_total
)


# Combine all data frames into one for ggplot
df_all <- rbind(df_month, df_hour, df_total)

# Plot using ggplot2 with faceting
figure3 <-ggplot(df_all, aes(x = Value, y = Fit)) +
  geom_line(color = "blue", size = 1) +  # Use geom_line instead of geom_smooth
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), 
              fill = rgb(0, 0, 1, 0.2), alpha = 0.3) + 
  facet_wrap(~ Variable, scales = "free_x") +  # Facet by variable
  labs(x = "Value", y = "Predicted Foraging Minutes", 
       title = "") +
  theme_minimal() +  # Minimal theme
  theme(strip.text = element_text(size = 14),  # Larger facet labels
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12),
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

figure3

#save plot
ggsave("Figure3_Foraging_DPM_Vessel totals.tiff", figure3, dpi=400)



########################
## AIC table
#######################

#run model 3 as if carrying on from model 2, removing the vars that we know are signif
#remove month first

model3<-bam(Feeding_Minutes~s(Hour, bs="cc",k=24) +s(Vessel_Total, bs="tp", k=6)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model3)

#now remove hour and add month back in 

model4<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc",k=12) +s(Vessel_Total, bs="tp", k=6)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model4)

#now remove vessel total and add hour back in 

model5<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc",k=12) +s(Hour, bs="cc", k=24)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model5)

#now remove year and add vessel total back in 

model6<-bam(Feeding_Minutes~s(Month_Numeric, bs="cc",k=12) +s(Hour, bs="cc", k=24)
            +s(Vessel_Total, bs="cc", k=6), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model6)


#run AICs
AIC(model2, model3,model4, model5,model6)
#df   AIC 
#model2 26.38844 36732.19
#model3 19.55749 37007.70
#model4 16.72218 37303.47
#model5 22.85733 36914.54
#model6 24.21406 36804.77

##################################


#####################
### NonFeeding minutes
######################

#VIF testing

testmodel2<- lm(NonFeeding_Minutes~Month_Numeric + Hour + Vessel_Total + Speed_Fast + Speed_Med + Speed_Slow, data = hp)
summary(testmodel2)

#test vif
vif(testmodel2)

              #GVIF^(1/(2*Df))
#Month_Numeric 1.006169 
#Hour          1.008505
#Vessel_Total  1.590953
#Speed_Fast    1.279216 
#Speed_Med     1.340535
#Speed_Slow    1.019583 

## Use GVIF^(1/2) to make comparable to other vif estimates
## none above 5 so nothing to drop

### TESTS FOR AUTOCORRELATION

#run a base model - without autocorrelation

hpbambase2<- bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24) +s(Vessel_Total, bs="tp", k=6) +s(Speed_Fast, bs="tp", k=5) 
                +s(Speed_Med, bs="tp",k=6) + s(Speed_Slow, bs="tp", k=4) +as.factor(Year), data=hp, family=negbin(0.1))
#inspect correlation
par(mfrow=c(1,1), cex=1.1)


# default ACF function:
acf(resid(hpbambase2), main="acf(resid(hpbambase2))")

#Determine the value of lag 1, as indicated by the red dot in the picture below:
# we also ask to plot the ACF by specifying plot (FALSE by default):
r1hp <- start_value_rho(hpbambase2, plot=TRUE)

#this gives the r1 value
acf(resid(hpbambase2), plot=FALSE)$acf[2]
#r1 = 0.4076227

### above 0.2 so need to run more tests


#determine Ar.start value
hp2 <- start_event(hp, column="Hour", event=c("Year","Date"),label.event="Event")

HP_new_Ar2 <-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24) +s(Vessel_Total, bs="tp", k=6) +s(Speed_Fast, bs="tp", k=5) 
                +s(Speed_Med, bs="tp",k=6) + s(Speed_Slow, bs="tp", k=4) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

#To retrieve the corrected residuals of the model, one could use the function resid_gam.
par(mfrow=c(1,2), cex=1.1)

# normal residuals:
normal.res <- resid(HP_new_Ar2)
acf(normal.res, main="HP normal")

# corrected residuals:
corrected.res <- resid_gam(HP_new_Ar2)
acf(corrected.res,main="HP corrected")


##### MODELLING TIME!

#model3

model3<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24) +s(Vessel_Total, bs="tp", k=6) +s(Speed_Fast, bs="tp", k=5) 
            +s(Speed_Med, bs="tp",k=6) + s(Speed_Slow, bs="tp", k=4) +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model3)

#Parametric coefficients:
                      #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          1.07614    0.18216   5.908 3.56e-09 ***
#as.factor(Year)2021 -0.23649    0.05230  -4.522 6.19e-06 ***
#as.factor(Year)2022 -0.49609    0.08131  -6.102 1.08e-09 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                  #edf Ref.df      F p-value    
#s(Month_Numeric) 8.379 10.000 33.645  <2e-16 ***
#s(Hour)          9.330 22.000 17.746  <2e-16 ***
#s(Vessel_Total)  4.206  4.633 37.851  <2e-16 ***
#s(Speed_Fast)    1.000  1.000  0.806   0.369    
#s(Speed_Med)     1.001  1.001  0.228   0.633    
#s(Speed_Slow)    1.001  1.001  1.575   0.209    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.198   Deviance explained = 16.3%
#fREML = 9614.2  Scale est. = 0.32388   n = 12754


#model4 (remove speed_fast, speed_med and speed_slow)

model4<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24) +s(Vessel_Total, bs="tp", k=6)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model4)

#Parametric coefficients:
                      #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.97536    0.06428  15.173  < 2e-16 ***
#as.factor(Year)2021 -0.23563    0.05240  -4.497 6.95e-06 ***
#as.factor(Year)2022 -0.49533    0.08148  -6.079 1.24e-09 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                #edf Ref.df     F p-value    
#s(Month_Numeric) 8.386 10.000 33.45  <2e-16 ***
#s(Hour)          9.313 22.000 17.87  <2e-16 ***
#s(Vessel_Total)  4.216  4.635 85.47  <2e-16 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.198   Deviance explained = 16.3%
#fREML = 9609.3  Scale est. = 0.32526   n = 12754

#run AICs between 2 models
AIC(model3,model4)

#df      AIC
#model3 29.65473 64736.93
#model4 26.63895 64730.59


#AIC for model4 marginally smaller and contains non-sig variables = better

##difference in explained deviance
1-(model3$dev/model3$null)
# 0.1628567

##difference in explained deviance
1-(model4$dev/model4$null)
# 0.1629228


#perform anova to compare models (start with model with least variables)
anova(model4, model3, test="Chisq")

#Resid. Df Resid. Dev       Df Deviance  Pr(>Chi)    
#1     12724     3954.9                         
#2     12721     3955.3 3.0309 -0.31239


### accept null hypothesis that all models are the same 

####### ALL VARIABLES TOGETHER IN FACET ########

# Refit the model after ensuring the factors are correctly defined
model4<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc", k=12) + s(Hour, bs="cc",k=24) +s(Vessel_Total, bs="tp", k=6),
             data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

# Create a data frame to store predictions for all variables

# Generate predictions for Month
new_data_month <- data.frame(
  Month_Numeric = seq(min(hp$Month_Numeric), max(hp$Month_Numeric), length.out = 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Total = rep(mean(hp$Vessel_Total), 300))

pred_month <- predict(model4, newdata = new_data_month, type = "response", se.fit = TRUE)
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
  Vessel_Total = rep(mean(hp$Vessel_Other), 300)
)

pred_hour <- predict(model4, newdata = new_data_hour, type = "response", se.fit = TRUE)
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

# Generate predictions for Vessel_Total
new_data_total <- data.frame(
  Month_Numeric = rep(mean(hp$Month_Numeric), 300),
  Hour = rep(mean(hp$Hour), 300),
  Vessel_Total = seq(min(hp$Vessel_Total), max(hp$Vessel_Total), length.out = 300)
)

pred_total <- predict(model4, newdata = new_data_total, type = "response", se.fit = TRUE)
fit_total <- pred_total$fit
se_total <- pred_total$se.fit
lower_ci_total <- fit_total - 1.96 * se_total
upper_ci_total <- fit_total + 1.96 * se_total

# Store predictions for Vessel_total
df_total <- data.frame(
  Variable = "Vessel_total",
  Value = new_data_total$Vessel_Total,
  Fit = fit_total,
  Lower_CI = lower_ci_total,
  Upper_CI = upper_ci_total
)


# Combine all data frames into one for ggplot
df_all <- rbind(df_month, df_hour, df_total)

# Plot using ggplot2 with faceting
figure4 <-ggplot(df_all, aes(x = Value, y = Fit)) +
  geom_line(color = "green", size = 1) +  # Use geom_line instead of geom_smooth
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), 
              fill = "#90EE90", alpha = 0.3) + 
  facet_wrap(~ Variable, scales = "free_x") +  # Facet by variable
  labs(x = "Value", y = "Predicted Foraging Minutes", 
       title = "") +
  theme_minimal() +  # Minimal theme
  theme(strip.text = element_text(size = 14),  # Larger facet labels
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12),
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

figure4

#save plot
ggsave("Figure4_Non-Foraging_DPM_Vessel totals.tiff", figure4, dpi=400)


########################
## AIC table
#######################

#run model 5 as if carrying on from model 4, removing the vars that we know are signif
#remove month first

model5<-bam(NonFeeding_Minutes~s(Hour, bs="cc",k=24) +s(Vessel_Total, bs="tp", k=6)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model5)

#now remove hour and add month back in 

model6<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc",k=12) +s(Vessel_Total, bs="tp", k=6)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model6)

#now remove vessel total and add hour back in 

model7<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc",k=12) +s(Hour, bs="cc",k=24)
            +as.factor(Year), data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model7)

#now remove year and add vessel total back in 

model8<-bam(NonFeeding_Minutes~s(Month_Numeric, bs="cc",k=12) +s(Hour, bs="cc",k=24)+s(Vessel_Total, bs="tp", k=6),
            data=hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model8)


#run AICs
AIC(model4, model5,model6, model7,model8)
#df   AIC 
#model4 26.63895 64730.59
#model5 18.63329 64933.60
#model6 17.07446 64947.92
#model7 22.97921 64920.54
#model8 24.60549 64752.57

##################################
## END
#################################


