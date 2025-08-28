# Foraging vs non foraging C-POD GAM models for visual vessels
# 18 Oct 2024
# Author: Dr. Chloe V. Robinson

RStudio.Version()
R.Version()

library(mgcv) #bam() function for models
library(car) #vif() function

#read in data

hp <- read.csv(("CPOD_vessels.csv"), stringsAsFactors =TRUE, na.strings = c(""))
attach(hp)
str(hp)

############################################
### Feeding trains with vessel categories
############################################

#VIF testing
#run without activity as this confounds other variables

testmodel1<- lm(Feeding_Trains~Month_Numeric + Hour + Type_FV + Type_F + Type_MVS + Type_ST
                + Type_T + Type_MVL + Type_SB, data = hp)
summary(testmodel1)

#test vif
vif(testmodel1)

                       #GVIF^(1/(2*Df))
#Month_Numeric        1.029504 
#Hour                 1.053816 
#Type_FV              1.017288
#Type_F               1.008670  
#Type_MVS             1.040006 
#Type_ST              1.049358
#Type_T               1.054011
#Type_MVL             1.022999 
#Type_SB              1.024544 

## Use GVIF^(1/2) to make comparable to other vif estimates
## none above 5 so nothing to drop

### TESTS FOR AUTOCORRELATION

#run a base model - without autocorrelation


hpbambase2<- bam(Feeding_Trains~Month_Numeric + Hour + Type_FV + Type_F + Type_MVS + Type_ST
                 + Type_T + Type_MVL + Type_SB, data = hp, family=negbin(0.1))
#inspect correlation
par(mfrow=c(1,1), cex=1.1)


# default ACF function:
acf(resid(hpbambase2), main="acf(resid(hpbambase2))")

#Determine the value of lag 1, as indicated by the red dot in the picture below:
# we also ask to plot the ACF by specifying plot (FALSE by default):
r1hp <- start_value_rho(hpbambase2, plot=TRUE)

#this gives the r1 value
acf(resid(hpbambase2), plot=FALSE)$acf[2]
#r1 = 0.1352938

## below 0.2 so no need to adjust



##### MODELLING TIME!

#model1

#had to remove  Type_F and Type_ST as too few values to support a K=1

model1<-bam(Feeding_Trains~s(Month_Numeric, bs="cc", k=9) + s(Hour, bs="cc",k=10) +s(Type_FV, bs="tp", k=5) 
             + s(Type_MVS, bs="tp",k=4) + s(Type_T, bs="tp", k=4)
            + s(Type_MVL, bs="tp", k=3) + s(Type_SB, bs="tp", k=3) +as.factor(Year), data = hp, family=negbin(0.1), discrete = TRUE)


summary(model1)

#Parametric coefficients:
                #Estimate Std. Error t value Pr(>|t|)   
#(Intercept)          -45.3622  8706.3156  -0.005  0.99584   
#as.factor(Year)2017   -2.3732     0.7410  -3.203  0.00138 **
#as.factor(Year)2018   -2.1142     0.7777  -2.718  0.00660 **
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                  #edf Ref.df     F  p-value    
#s(Month_Numeric) 3.752e+00   7.00 3.741 6.01e-06 ***
#s(Hour)          6.445e-06   8.00 0.000    0.622    
#s(Type_FV)       1.000e+00   1.00 1.466    0.226    
#s(Type_MVS)      1.000e+00   1.00 1.345    0.246    
#s(Type_T)        1.190e+00   1.12 0.000    1.000    
#s(Type_MVL)      1.000e+00   1.00 1.329    0.249    
#s(Type_SB)       1.000e+00   1.00 0.000    0.998    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.0308   Deviance explained = 35.7%
#fREML = 4093.4  Scale est. = 2.5661    n = 2911


#model2 (remove non-sig values)

model2<-bam(Feeding_Trains~s(Month_Numeric, bs="cc", k=9) +as.factor(Year), data = hp, family=negbin(0.1), discrete = TRUE)
summary(model2)

#Parametric coefficients:
                        #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          -1.5854     0.5616  -2.823  0.00479 **
#as.factor(Year)2017  -2.4681     0.7969  -3.097  0.00197 **
#as.factor(Year)2018  -2.0845     0.8665  -2.406  0.01620 * 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                #edf Ref.df     F  p-value    
#s(Month_Numeric) 3.273      7 2.658 0.000156 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.0308   Deviance explained =   30%
#fREML = 4645.7  Scale est. = 3.7422    n = 2911


#run AICs between 2 models
AIC(model1,model2)

#df      AIC
#model1 11.919181 667.5528
#model2  7.743863 686.5155


#AIC for model1 is lower but has non sig variables (best to use model2 as has fewer variables)

##difference in explained deviance
1-(model1$dev/model1$null)
#0.3568585

##difference in explained deviance
1-(model2$dev/model2$null)
#0.2995774


#perform anova to compare models (start with model with least variables)
anova(model2, model1, test="Chisq")

#Resid. Df Resid. Dev       Df Deviance  Pr(>Chi)    
#1     2903.4     333.98                           
#2    2898.7     306.67 4.7155   27.313  0.04968 *


### reject null hypothesis that both models are the same 

########################################



#####################################
### Feeding trains with vessel activity
######################################

#VIF testing
#run without types as this confounds other variables

testmodel2<- lm(Feeding_Trains~Month_Numeric + Hour + Activity_Fishing + Activity_Stopped + Activity_Transit +
               Activity_Work + Year, data = hp)
summary(testmodel2)

#test vif
vif(testmodel2)

#GVIF^(1/(2*Df))
#Month_Numeric        1.151304  
#Hour                 1.065072 
#Activity_Fishing     1.013468
#Activity_Stopped     1.033940  
#Activity_Transit     1.020262 
#Activity_Work        1.022916
#Year                1.210122 
 

## Use GVIF^(1/2) to make comparable to other vif estimates
## none above 5 so nothing to drop

### TESTS FOR AUTOCORRELATION

#run a base model - without autocorrelation


hpbambase2<- bam(Feeding_Trains~Month_Numeric + Hour + Activity_Fishing + Activity_Stopped + Activity_Transit +
                   Activity_Work + Year, data = hp, family=negbin(0.1))
#inspect correlation
par(mfrow=c(1,1), cex=1.1)


# default ACF function:
acf(resid(hpbambase2), main="acf(resid(hpbambase2))")

#Determine the value of lag 1, as indicated by the red dot in the picture below:
# we also ask to plot the ACF by specifying plot (FALSE by default):
r1hp <- start_value_rho(hpbambase2, plot=TRUE)

#this gives the r1 value
acf(resid(hpbambase2), plot=FALSE)$acf[2]
#r1 = 0.1180417

## below 0.2 so no need to adjust



##### MODELLING TIME!

#model1

#had to remove Activity_Stopped as too few values to support a K=1

model1<-bam(Feeding_Trains~s(Month_Numeric, bs="cc", k=9) + s(Hour, bs="cc",k=10) +s(Activity_Transit, bs="tp", k=6) 
            + s(Activity_Fishing, bs="tp",k=4) + s(Activity_Work, bs="tp", k=3) +as.factor(Year), data = hp, family=negbin(0.1), discrete = TRUE)


summary(model1)

#Parametric coefficients:
                      #Estimate Std. Error t value Pr(>|t|)   
#(Intercept)         -3.137e+01  1.099e+06   0.000  0.99998   
#as.factor(Year)2017 -2.343e+00  7.297e-01  -3.211  0.00134 **
#as.factor(Year)2018 -2.120e+00  7.568e-01  -2.801  0.00513 **
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
#edf Ref.df     F  p-value    
#s(Month_Numeric)    3.707e+00 7.0000 3.616 8.87e-06 ***
#s(Hour)             4.708e-06 8.0000 0.000    0.517    
#s(Activity_Transit) 1.000e+00 1.0000 1.757    0.185    
#s(Activity_Fishing) 1.657e+00 1.9626 0.636    0.554    
#s(Activity_Work)    9.922e-01 0.9946 0.000    0.500    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.0282   Deviance explained = 33.4%
#fREML = 4137.7  Scale est. = 2.6351    n = 2911


#model2 (remove non-sig values)

model2<-bam(Feeding_Trains~s(Month_Numeric, bs="cc", k=9) +as.factor(Year), data = hp, family=negbin(0.1), discrete = TRUE)
summary(model2)

#Parametric coefficients:
                    #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          -1.5854     0.5616  -2.823  0.00479 **
#as.factor(Year)2017  -2.4681     0.7969  -3.097  0.00197 **
#as.factor(Year)2018  -2.0845     0.8665  -2.406  0.01620 * 
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                #edf Ref.df     F  p-value    
#s(Month_Numeric) 3.273      7 2.658 0.000156 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.0308   Deviance explained =   30%
#fREML = 4645.7  Scale est. = 3.7422    n = 2911


#run AICs between 2 models
AIC(model1,model2)

#df      AIC
#model1 11.781351 678.3142
#model2  7.743863 686.5155


#AIC for model1 is lower but has non sig variables (best to use model2 as has fewer variables)

##difference in explained deviance
1-(model1$dev/model1$null)
#0.3337117

##difference in explained deviance
1-(model2$dev/model2$null)
#0.2995774


#perform anova to compare models (start with model with least variables)
anova(model2, model1, test="Chisq")

#Resid. Df Resid. Dev       Df Deviance  Pr(>Chi)    
#1     2903.4     333.98                         
#2    2898.9     317.71 4.4996   16.276   0.2358


### accept null hypothesis that all models are the same 

########################################

####### ALL VARIABLES TOGETHER IN FACET ########

# Refit the model after ensuring the factors are correctly defined (remove year for plotting smooth only)
model3<-bam(Feeding_Trains~s(Month_Numeric, bs="cc", k=9), data = hp, family=negbin(0.1), discrete = TRUE)

# Create a data frame to store predictions for all variables

# Generate predictions for Month
new_data_month <- data.frame(
  Month_Numeric = seq(min(hp$Month_Numeric), max(hp$Month_Numeric), length.out = 300))

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


dev.off()

# Plot using ggplot2 with faceting
figureS2 <-ggplot(df_month, aes(x = Value, y = Fit)) +
  geom_line(color = "pink", size = 1) +  # Use geom_line instead of geom_smooth
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), 
              fill = "lightgrey", alpha = 0.3) + 
  facet_wrap(~ Variable, scales = "free_x") +  # Facet by variable
  labs(x = "Value", y = "Predicted Foraging Trains", 
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

figureS2

#save plot
ggsave("Figures2_Foraging_train_mont_CPOD.tiff", figureS2, dpi=400)



################################################
### Non-Foraging trains with vessel categories
################################################

#VIF testing
#run without activity as this confounds other variables

testmodel1<- lm(NonFeeding_Trains~Month_Numeric + Hour + Type_FV + Type_F + Type_MVS + Type_ST
                + Type_T + Type_MVL + Type_SB, data = hp)
summary(testmodel1)

#test vif
vif(testmodel1)

#GVIF^(1/(2*Df))
#Month_Numeric        1.029504 
#Hour                 1.053816 
#Type_FV              1.017288
#Type_F               1.008670  
#Type_MVS             1.040006 
#Type_ST              1.049358
#Type_T               1.054011
#Type_MVL             1.022999 
#Type_SB              1.024544 

## Use GVIF^(1/2) to make comparable to other vif estimates
## none above 5 so nothing to drop

### TESTS FOR AUTOCORRELATION

#run a base model - without autocorrelation


hpbambase2<- bam(NonFeeding_Trains~Month_Numeric + Hour + Type_FV + Type_F + Type_MVS + Type_ST
                 + Type_T + Type_MVL + Type_SB, data = hp, family=negbin(0.1))
#inspect correlation
par(mfrow=c(1,1), cex=1.1)


# default ACF function:
acf(resid(hpbambase2), main="acf(resid(hpbambase2))")

#Determine the value of lag 1, as indicated by the red dot in the picture below:
# we also ask to plot the ACF by specifying plot (FALSE by default):
r1hp <- start_value_rho(hpbambase2, plot=TRUE)

#this gives the r1 value
acf(resid(hpbambase2), plot=FALSE)$acf[2]
#r1 = 0.3946055

## below 0.2 so need to adjust

#determine Ar.start value
hp <- start_event(hp, column="Hour", event=c("Year","Date"),label.event="Event")

HP_new_Ar2 <-bam(NonFeeding_Trains~s(Month_Numeric, bs="cc", k=9) + s(Hour, bs="cc",k=10) +s(Type_FV, bs="tp", k=5) 
                 + s(Type_MVS, bs="tp",k=4) + s(Type_T, bs="tp", k=4)
                 + s(Type_MVL, bs="tp", k=3) + s(Type_SB, bs="tp", k=3) +as.factor(Year), data = hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

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

#model1

#had to remove  Type_F and Type_ST as too few values to support a K=1

model1<-bam(NonFeeding_Trains~s(Month_Numeric, bs="cc", k=9) + s(Hour, bs="cc",k=10) +s(Type_FV, bs="tp", k=5) 
            + s(Type_MVS, bs="tp",k=4) + s(Type_T, bs="tp", k=4)
            + s(Type_MVL, bs="tp", k=3) + s(Type_SB, bs="tp", k=3) +as.factor(Year), data = hp, family=negbin(0.1),discrete = TRUE, rho=r1hp, AR.start = hp$start.event)


summary(model1)

#Parametric coefficients:
                      #Estimate Std. Error t value Pr(>|t|)   
#(Intercept)          -9.9399    31.8178  -0.312 0.754757    
#as.factor(Year)2017  -3.9919     0.8219  -4.857 1.26e-06 ***
#as.factor(Year)2018  -2.6841     0.7305  -3.674 0.000243 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                  #edf Ref.df     F  p-value    
#s(Month_Numeric) 4.973e+00  7.000 4.485 2.15e-06 ***
#s(Hour)          5.372e-06  8.000 0.000    0.692    
#s(Type_FV)       1.000e+00  1.000 1.415    0.234    
#s(Type_MVS)      1.000e+00  1.000 0.023    0.879    
#s(Type_T)        1.090e+00  1.172 0.299    0.512    
#s(Type_MVL)      1.000e+00  1.000 0.172    0.678    
#s(Type_SB)       1.000e+00  1.000 0.052    0.820    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.0612   Deviance explained = 41.5%
#fREML = 3654.9  Scale est. = 2.0133    n = 2911


#model2 (remove non-sig values)

model2<-bam(NonFeeding_Trains~s(Month_Numeric, bs="cc", k=9) +as.factor(Year), data = hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model2)

#Parametric coefficients:
                  #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          -0.3726     0.4579  -0.814 0.415965    
#as.factor(Year)2017  -4.0931     0.8172  -5.009 5.81e-07 ***
#as.factor(Year)2018  -2.6827     0.7270  -3.690 0.000228 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                #edf Ref.df     F  p-value    
#s(Month_Numeric) 4.885      7 5.188 2.09e-08 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.0547   Deviance explained = 39.7%
#fREML = 3718.7  Scale est. = 2.0963    n = 2911

#run AICs between 2 models
AIC(model1,model2)

#df      AIC
#model1 15.056434 992.2518
#model2  9.782201 995.1881


#AIC for model1 is lower but has more non-sig variables

##difference in explained deviance
1-(model1$dev/model1$null)
#0.4152692

##difference in explained deviance
1-(model2$dev/model2$null)
#0.3974957


#perform anova to compare models (start with model with least variables)
anova(model2, model1, test="Chisq")

    #Resid. Df Resid. Dev       Df Deviance  Pr(>Chi)    
#1    2901.3     457.12                         
#2    2895.9     443.64 5.3712   13.485   0.2821 


### accept null hypothesis that all models are the same 

########################################


####### ALL VARIABLES TOGETHER IN FACET ########

# Refit the model after ensuring the factors are correctly defined (remove year for plotting smooth only)
model3<-bam(NonFeeding_Trains~s(Month_Numeric, bs="cc", k=9), data = hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

# Create a data frame to store predictions for all variables

# Generate predictions for Month
new_data_month <- data.frame(
  Month_Numeric = seq(min(hp$Month_Numeric), max(hp$Month_Numeric), length.out = 300))

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


dev.off()

# Plot using ggplot2 with faceting
figureS3 <-ggplot(df_month, aes(x = Value, y = Fit)) +
  geom_line(color = "purple", size = 1) +  # Use geom_line instead of geom_smooth
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), 
              fill = "lightgrey", alpha = 0.3) + 
  facet_wrap(~ Variable, scales = "free_x") +  # Facet by variable
  labs(x = "Value", y = "Predicted Foraging Trains", 
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

figureS3

#save plot
ggsave("Figures3_Non-Foraging_train_month_CPOD.tiff", figureS3, dpi=400)


##########################################
### NonFeeding trains with vessel activity
##########################################

#VIF testing
#run without types as this confounds other variables

testmodel2<- lm(NonFeeding_Trains~Month_Numeric + Hour + Activity_Fishing + Activity_Stopped + Activity_Transit +
                  Activity_Work + Year, data = hp)
summary(testmodel2)

#test vif
vif(testmodel2)

#GVIF^(1/(2*Df))
#Month_Numeric        1.151304  
#Hour                 1.065072 
#Activity_Fishing     1.013468
#Activity_Stopped     1.033940  
#Activity_Transit     1.020262 
#Activity_Work        1.022916
#Year                1.210122 


## Use GVIF^(1/2) to make comparable to other vif estimates
## none above 5 so nothing to drop

### TESTS FOR AUTOCORRELATION

#run a base model - without autocorrelation


hpbambase2<- bam(NonFeeding_Trains~Month_Numeric + Hour + Activity_Fishing + Activity_Stopped + Activity_Transit +
                   Activity_Work + Year, data = hp, family=negbin(0.1))
#inspect correlation
par(mfrow=c(1,1), cex=1.1)


# default ACF function:
acf(resid(hpbambase2), main="acf(resid(hpbambase2))")

#Determine the value of lag 1, as indicated by the red dot in the picture below:
# we also ask to plot the ACF by specifying plot (FALSE by default):
r1hp <- start_value_rho(hpbambase2, plot=TRUE)

#this gives the r1 value
acf(resid(hpbambase2), plot=FALSE)$acf[2]
#r1 = 0.4249391

## above 0.2 so need to adjust

#determine Ar.start value
hp <- start_event(hp, column="Hour", event=c("Year","Date"),label.event="Event")

HP_new_Ar2 <-bam(NonFeeding_Trains~s(Month_Numeric, bs="cc", k=9) + s(Hour, bs="cc",k=10) +s(Activity_Transit, bs="tp", k=6) 
                 + s(Activity_Fishing, bs="tp",k=4) + s(Activity_Work, bs="tp", k=3) +as.factor(Year), data = hp, family=negbin(0.1), discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

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

#model1

#had to remove Activity_Stopped as too few values to support a K=1

model1<-bam(NonFeeding_Trains~s(Month_Numeric, bs="cc", k=9) + s(Hour, bs="cc",k=10) +s(Activity_Transit, bs="tp", k=6) 
            + s(Activity_Fishing, bs="tp",k=4) + s(Activity_Work, bs="tp", k=3) +as.factor(Year), data = hp, family=negbin(0.1),discrete = TRUE, rho=r1hp, AR.start = hp$start.event)

summary(model1)

#Parametric coefficients:
                  #Estimate Std. Error t value Pr(>|t|)
#(Intercept)          -0.6290     2.1122  -0.298  0.76588    
#as.factor(Year)2017  -4.0754     0.8269  -4.928 8.76e-07 ***
#as.factor(Year)2018  -2.7177     0.7160  -3.796  0.00015 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                    #edf Ref.df     F  p-value    
#s(Month_Numeric)    4.964e+00  7.000 4.633 1.21e-06 ***
#s(Hour)             8.222e-06  8.000 0.000   0.0285 *  
#s(Activity_Transit) 1.542e+00  1.884 0.769   0.4237    
#s(Activity_Fishing) 1.415e+00  1.686 0.550   0.5156    
#s(Activity_Work)    1.000e+00  1.000 0.066   0.7972    
#---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.0548   Deviance explained = 41.4%
#fREML = 3571.8  Scale est. = 1.9237    n = 2911


#model2 (remove non-sig values)

model2<-bam(NonFeeding_Trains~s(Month_Numeric, bs="cc", k=9) + s(Hour, bs="cc",k=10) +as.factor(Year), data = hp, family=negbin(0.1),discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model2)

#Parametric coefficients:
                    #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          -0.3025     0.5664  -0.534 0.593330    
#as.factor(Year)2017  -4.0142     0.8381  -4.790 1.75e-06 ***
#as.factor(Year)2018  -2.5735     0.7352  -3.500 0.000472 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                #edf Ref.df    F  p-value    
#s(Month_Numeric) 4.659e+00      7 4.27 3.16e-06 ***
#s(Hour)          1.103e-05      8 0.00    0.199    
#---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.052   Deviance explained = 38.7%
#fREML =   3790  Scale est. = 2.2547    n = 2911


#model3 (remove non-sig values)

model3<-bam(NonFeeding_Trains~s(Month_Numeric, bs="cc", k=9)
             +as.factor(Year), data = hp, family=negbin(0.1),discrete = TRUE, rho=r1hp, AR.start = hp$start.event)
summary(model3)

#Parametric coefficients:
                    #Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          -0.4184     0.4723  -0.886 0.375858    
#as.factor(Year)2017  -4.0811     0.8458  -4.825 1.47e-06 ***
#as.factor(Year)2018  -2.6442     0.7481  -3.534 0.000415 ***
  #---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
                  #edf Ref.df    F  p-value    
#s(Month_Numeric) 4.623      7 4.66 9.47e-07 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#R-sq.(adj) =  0.0532   Deviance explained = 39.2%
#fREML = 3785.5  Scale est. = 2.2474    n = 2911


#run AICs between 2 models
AIC(model1,model2, model3)

#df      AIC
#model1 14.439020 937.4797
#model2  9.586046 948.7626
#model3  9.542950 944.2103


#AIC for model3 is lowest but has non-sig variables 


##difference in explained deviance
1-(model1$dev/model1$null)
#0.4142187

##difference in explained deviance
1-(model2$dev/model2$null)
#0.3865545

##difference in explained deviance
1-(model3$dev/model3$null)
#0.3924411


#perform anova to compare models (start with model with least variables)
anova(model3, model2, model1, test="Chisq")

      #Resid. Df Resid. Dev       Df Deviance  Pr(>Chi)    
#1    2901.5     460.96                             
#2    2901.5     465.42 0.050427  -4.4662           
#3    2896.0     444.43 5.444316  20.9888  0.06845 .


### accept null hypothesis that all models are the same 

########################################
## END
#######################################
