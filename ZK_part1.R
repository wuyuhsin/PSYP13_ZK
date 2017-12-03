####################################################################
#                                                                  #
#                         ZK. Assignment 1                         #
#                                                                  #
####################################################################


# loading packages
library(psych) # for describe
library(lm.beta) # for lm.beta
library(dplyr)
library(gsheet)
library(car) # for scatter3d, residualPlots, vif, pairs.panels, ncvTest
library(ggplot2) # for ggplot
library(rgl) # for scatter3d
library(lsr) # for standardCoefs


####################################################################
#                Load dataset and explore data                     #
####################################################################

# import data set
data_sample_1 <-  read.csv("https://raw.githubusercontent.com/kekecsz/PSYP13_Data_analysis_class/master/home_sample_1.csv")

# descriptive statistics
View(data_sample_1)
summary(data_sample_1)
describe(data_sample_1)

# ckecking the coding errors
which(data_sample_1$sex==3) # error that was typed 3 in the sex variable (ID_15)
which(data_sample_1$mindfulness<0) # value below 0 in the mindfulness variable (ID_24)

# delete the case ID_15 and ID_24 having coding error
data_sample_1_delete.coding.error <- data_sample_1[-c(15,24),]
summary(data_sample_1_delete.coding.error)
describe(data_sample_1_delete.coding.error)

# histograms
par(mfcol=c(2,2))
hist(data_sample_1_delete.coding.error$pain, main = "", xlab = "pain")
hist(data_sample_1_delete.coding.error$age, main = "", xlab = "age")
hist(data_sample_1_delete.coding.error$STAI_trait, main = "", xlab = "STAI_trait")
hist(data_sample_1_delete.coding.error$pain_cat, main = "", xlab = "pain_cat")
hist(data_sample_1_delete.coding.error$cortisol_serum, main = "", xlab = "cortisol_serum")
hist(data_sample_1_delete.coding.error$cortisol_saliva, main = "", xlab = "cortisol_saliva")
hist(data_sample_1_delete.coding.error$mindfulness, main = "", xlab = "mindfulness")
hist(data_sample_1_delete.coding.error$weight, main = "", xlab = "weight")

# scatterplot matrix
par(mfcol=c(1,1))
pairs(data_sample_1_delete.coding.error[,-c(1,3)], pch = ".", cex = 1.5)

# correlation matrix
round(cor(data_sample_1_delete.coding.error[,-c(1,3)]), 4) 


####################################################################
#                    Build regression models                       #
####################################################################

mod1 <- lm( pain ~ sex + age, 
            data = data_sample_1_delete.coding.error)

mod2 <- lm( pain ~ sex + age + STAI_trait + pain_cat + mindfulness + cortisol_serum + cortisol_saliva, 
            data = data_sample_1_delete.coding.error)



####################################################################
#                  Data and model diagnostics                      #
####################################################################

### check the influential outliers

# Cook’s distance
cook <- cooks.distance ( model = mod2 )
which(cook > 1)  # cut off >1
plot(mod2, which = 4)

# residuals vs leverage
plot(mod2, which = 5)

# Mahalanobois distance
lev <- hat(model.matrix(mod2))
data_sample_1_delete.coding.error[lev > .045, ] 

N <- nrow(data_sample_1_delete.coding.error)
mahad <- (N-1) * (lev-1/N)
tail(sort(mahad), 7)
order(mahad, decreasing = TRUE)[c(7,6,5,4,3,2,1)] # cut off >24.32



### check assumptions

# 1) normality assumption
# QQ plot
plot(mod2, which = 2)
# skew and kurtosis
describe(residuals(mod2))
# histogram
hist(residuals(mod2), breaks = 20)

# 2) linearity assumption
# predicted values against actual values
pred <- predict( object = mod2 )
plot( x = pred, y = data_sample_1_delete.coding.error$pain, 
      xlab = "Fitted Values", ylab = "Observed Values")

# predicted values against residuals
plot(mod2, which = 1)

# residual plot for each predictor from the car package, returning the result of linearity tests
residualPlots(mod2)


# 3) homoscedasticty assumption (homogeneity of variance)
plot(mod2, which = 3)
ncvTest(mod2)

# 4) multicollinearity (VIF above 5), or a VIF threshold of 3 is recommended
vif(mod2)
pairs.panels(data_sample_1_delete.coding.error[,c("pain", "sex", "age","STAI_trait","pain_cat","mindfulness","cortisol_serum","cortisol_saliva")], 
             col = "red", lm = T)



####################################################################
#                    Hierarchical regression                       #
####################################################################

# Quantify the amount of information gained about the variability of the outcome
# by adding pedictors in block certain predictors

mod1 <- lm( pain ~ sex + age, 
            data = data_sample_1_delete.coding.error)

mod2 <- lm( pain ~ sex + age + STAI_trait + pain_cat + mindfulness + cortisol_serum + cortisol_saliva, 
            data = data_sample_1_delete.coding.error)

summary(mod1)
summary(mod2)


# 95% confidence intervals for the coefficient
confint(object = mod1)
confint(object = mod2)

# standardised regression coefficient
standardCoefs(mod1)
standardCoefs(mod2)



####################################################################
#                       Model comparisons                          #
####################################################################

# compare models with the anova function
anova(mod1, mod2)

# compare with AIC
AIC(mod1)
AIC(mod2)





