####################################################################
#                                                                  #
#                         ZK. Assignment 2                         #
#                                                                  #
####################################################################


# load packages
library(psych) # for describe
library(car) # for residualPlots, vif, pairs.panels, ncvTest
library(ggplot2) # for ggplot
library(cAIC4) # for cAIC
library(r2glmm) # for r2beta
library(influence.ME) # for influence
library(lattice) # for qqmath
library(reshape2) # for melt function
library(lsr)


# coustom functions
# function to extract standardized beta coefficients from mer models, such as one produced by lmer
# source: https://stackoverflow.com/questions/25142901/standardized-coefficients-for-lmer-model

stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}

####################################################################
#                          Load dataset                            #
####################################################################

# load data
data_sample_1 <-  read.csv("https://raw.githubusercontent.com/kekecsz/PSYP13_Data_analysis_class/master/home_sample_1.csv")
summary(data_sample_1)

# delete the coding error object
data_sample_1_delete.coding.error <- data_sample_1[-c(15,24),]
summary(data_sample_1_delete.coding.error)


####################################################################
#              New model suggested by the oppoent                  #
####################################################################

# the oppoent used the following variables as predictors in the initial model : 
# age, sex, weight, STAI, pain catastrophizing, mindfulness, serum cortisol. 

mod_oppoent <- lm(pain ~ sex + age + weight + STAI_trait + pain_cat + mindfulness + cortisol_serum, 
                  data = data_sample_1_delete.coding.error)

summary(mod_oppoent)


####################################################################
#             Re-run the data and model diagnostics                #
####################################################################

## checking for influential outliers

# Cook's distance
cook <- cooks.distance ( model = mod_oppoent )
which(cook > 1) 
plot(mod_oppoent, which = 4) 

# residuals vs leverage
plot(mod_oppoent, which = 5) 

# Mahalanobois distance
lev <- hat(model.matrix(mod_oppoent))
data_sample_1_delete.coding.error[lev > .045, ] 

N <- nrow(data_sample_1_delete.coding.error)
mahad <- (N-1) * (lev-1/N)
tail(sort(mahad), 7)
order(mahad, decreasing = TRUE)[c(7,6,5,4,3,2,1)] # cut off >24.32 



## checking assumptions

# 1) normality assumption
# QQ plot
plot(mod_oppoent, which = 2)
# skew and kurtosis
describe(residuals(mod_oppoent))
# histogram
hist(residuals(mod_oppoent), breaks = 20)


# 2) linearity assumption
pred <- predict( object = mod_oppoent )
plot( x = pred, y = data_sample_1_delete.coding.error$pain, 
      xlab = "Fitted Values", ylab = "Observed Values")
# predicted values against residuals
plot(mod_oppoent, which = 1)
# residual plot for each predictor from the car package, returning the result of linearity tests
residualPlots(mod_oppoent)


# 3) homoscedasticty assumption (homogeneity of variance)
plot(mod_oppoent, which = 3)
ncvTest(mod_oppoent)


# 4) multicollinearity (VIF above 5), or a VIF threshold of 3 
vif(mod_oppoent)
pairs.panels(data_sample_1_delete.coding.error[,c("pain", "sex", "age", "weight", "STAI_trait", "pain_cat", "mindfulness", "cortisol_serum")], col = "red", lm = T)



####################################################################
#                       Backward regression                        #
####################################################################

# backward regression to identify the variables with the highest unique predictive value
mod_oppoent_back <-  step(mod_oppoent, direction = "backward")

# model with predictors that were retained in the backward regression: age, pain_cat, mindfulness, cortisol_saliva 
backward_model <- lm(pain ~ age + pain_cat + mindfulness + cortisol_serum, 
                              data = data_sample_1_delete.coding.error)

# compute CI 
confint(backward_model) 

# compute the standardised regression coefficient
standardCoefs(backward_model)

# model in assignment 1
theory_based_model <- lm(pain ~ sex + age + STAI_trait + pain_cat + mindfulness + cortisol_serum + cortisol_saliva, 
                         data = data_sample_1_delete.coding.error)



# compare the initial model and the backward model
summary(mod_oppoent)
summary(backward_model)

AIC(mod_oppoent)
AIC(backward_model)

anova(backward_model, mod_oppoent)


# compare the backward model and the theory based model
summary(backward_model)
summary(theory_based_model)

AIC(backward_model) 
AIC(theory_based_model) 

anova(backward_model, theory_based_model) 


####################################################################
#         Put the models to the test on some new data              #
####################################################################

# load new data from another 160 participants 
data_sample_2 <-read.csv("https://raw.githubusercontent.com/kekecsz/PSYP13_Data_analysis_class/master/home_sample_2.csv")
View(data_sample_2)
summary(data_sample_2)

# make predictions on pain using the backward model which were trained on data file 1
pred_test_backward <- predict(backward_model, data_sample_2)

# make predictions on pain using the theory-based model which were trained on data file 1
pred_test_theory <- predict(theory_based_model, data_sample_2)

# calculate the sum of squared residuals 
RSS_test_backward <-  sum((data_sample_2[,"pain"] - pred_test_backward)^2)
RSS_test_theory <-  sum((data_sample_2[,"pain"] - pred_test_theory)^2)
RSS_test_backward 
RSS_test_theory 







