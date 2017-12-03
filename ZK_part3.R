####################################################################
#                                                                  #
#                         ZK. Assignment 3                         #
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
library(lmerTest)

                                                        
# Custom functions:                        
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
#                           Load dataset                           #
####################################################################

# load dataset
data_sample_3 <- read.csv("https://raw.githubusercontent.com/kekecsz/PSYP13_Data_analysis_class/master/home_sample_3.csv")
View(data_sample_3)

# asign ID and sex as factors
data_sample_3$ID <- factor(data_wound$ID)
data_wound$sex <- factor(data_wound$sex)

# variables
names(data_sample_3)

# designate the repeated varibales
repeated_variables <- c("pain1",	"pain2", "pain3",	"pain4")


####################################################################
#                            Explore data                          #
####################################################################

# descriptives
describe(data_sample_3)

# histograms
par(mfcol=c(2,2))
hist(data_sample_3$pain1)
hist(data_sample_3$pain2)
hist(data_sample_3$pain3)
hist(data_sample_3$pain4)

# correlation of repeated variables
cor(data_sample_3[,repeated_variables])


####################################################################
#        Convert the dataframe from wide to long format            #
####################################################################

# ID vars should be all non-repeated variables
data_sample_3_long <-  melt(data_sample_3, 
                            measure.vars = repeated_variables, 
                            variable.name = "day", 
                            value.name = "pain")

# order data frame by participant ID
data_sample_3_long <- data_sample_3_long[order(data_sample_3_long[,"ID"]),]

# change the time variable to a numerical vector
data_sample_3_long$day <- as.numeric(data_sample_3_long$day)

# data in long format
data_sample_3_long
View(data_sample_3_long)


####################################################################
#                  Fit linear mixed effect models                  #
####################################################################

# random intercept model
mod_rep_int <-  lmer(pain ~ age + sex + weight + STAI_trait + pain_cat + mindfulmess + cortisol_serum + day + (1|ID), 
                     data = data_sample_3_long)
summary(mod_rep_int)

# random slope and intercept model
mod_rep_slope <- lmer(pain ~ age + sex + weight + STAI_trait + pain_cat + mindfulness + cortisol_serum + day + (day|ID), 
                      data = data_sample_3_long)
summary(mod_rep_slope)


####################################################################
#                       Model comparisons                          #
####################################################################

## plot the pain ratings (y axis) over time (x axis) for each participant separately 
## in separate panels (or facets) in a single graph

# save the predictions of both models to variables
data_sample_3_long$pred_int <-  predict(mod_rep_int)
data_sample_3_long$pred_slope <-  predict(mod_rep_slope)

# plot the random intercept model
ggplot(data_sample_3_long, aes(y = pain, x = day, group = ID))+
  geom_point(size = 2)+
  geom_line(color = 'red', aes(y = pred_int, x = day))+
  facet_wrap( ~ ID, ncol = 5)

# plot the random slope and intercept model 
ggplot(data_sample_3_long, aes(y = pain, x = day, group = ID))+
  geom_point(size = 2)+
  geom_line(color = 'red', aes(y = pred_slope, x = day))+
  facet_wrap( ~ ID, ncol = 5)

## compare models with cAIC
cAIC(mod_rep_int)$caic 
cAIC(mod_rep_slope)$caic 

## compare models with anova
anova(mod_rep_int, mod_rep_slope) 



####################################################################
#           Add a quadratic term of time to the model              #
####################################################################

## random slope and intercept model including a quadratic term of time:
## (account for curved relationship between time and pain)
mod_rep_slope_quad <-  lmer(pain ~ age + sex + weight + STAI_trait + pain_cat + mindfulness + cortisol_serum + day + (day|ID) + I(day^2), 
                          data = data_sample_3_long)

## plot the results
# save prediction of the model to new variable
data_sample_3_long$pred_slope_quad <-  predict(mod_rep_slope_quad)

# plot the resulting regression line over the actual pain ratings of each participants
ggplot(data_sample_3_long, aes(y = pain, x = day, group = ID))+
  geom_point(size = 2)+
  geom_line(color='red', aes(y = pred_slope_quad, x = day))+
  facet_wrap( ~ ID, ncol = 5)

## compare models with cAIC
cAIC(mod_rep_slope)$caic 
cAIC(mod_rep_slope_quad)$caic 

## compare models with anova
anova(mod_rep_slope, mod_rep_slope_quad) 



####################################################################
#                 Final model: mod_rep_slope_quad                  #
####################################################################

summary(mod_rep_slope_quad)

# compute the marginal R2 (proportion of variance explained by the fixed factors alone)
r2beta(mod_rep_slope_quad, method = "nsj", data = data_sample_3_long)  # (method='nsj', the Nakagawa and Schielzeth approach is applied)

# compute the CI 
confint(mod_rep_slope_quad)

# compute the Beta (standardized coefficients)
stdCoef.merMod(mod_rep_slope_quad)



####################################################################
#                       Model diagnostics                          #
####################################################################

# 1) normality assumption
# QQ plot
qqmath(mod_rep_slope_quad, id=0.05) 


# 2) linearity assumption
# linearity of prediction and standardized residuals
plot(mod_rep_slope_quad)
# linearity of each predictor and the standardized residual
# visualize the linearity of connection of each predictor
predictors <- c("age", "sex", "weight", "STAI_trait", "pain_cat", "cortisol_serum", "day")

for(i in 1:length(predictors)){
  predictor_to_test = data_sample_3_long[,predictors[i]]
  print(
    ggplot(data.frame(x = predictor_to_test, pearson = residuals(mod_rep_slope_quad, type = "pearson")),
           aes(x = x, y = pearson)) +
      geom_point() +
      geom_smooth(method = 'loess') +
      theme_bw()
  )
}


# 3) homoscedasticty assumption (homogeneity of variance)
plot(mod_rep_slope_quad)


# 4) multicollinearity
pairs.panels(data_sample_3_long, col = "red", lm = T)







