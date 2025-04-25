#### Load necessary packages and data ####

# Loading required libraries
library(tidyverse)  # Data manipulation and visualization
library(brms)       # Bayesian regression models
library(emmeans)    # Estimated marginal means
library(lme4)       # Linear mixed-effects models
library(ggpubr)     # Publication-ready plots
library(MuMIn)      # Model selection and averaging
library(paletteer)  # Color palettes
library(car)        # Companion to applied regression
library(lmerTest)   # Test for mixed-effects models
library(performance) # Model performance checks
library(datawizard) # Data transformations (standardization)
library(DHARMa)     # Residual diagnostics for mixed models
library(glmmTMB)    # Generalized linear mixed-effects models

# Set working directory to where the data is located
setwd("C:/Users/alice/Documents/PostDoc/Fieldwork/Playback/Analysis")

# Load the dataset
data <- read.csv(file = "Pb_data_clean.csv")

# Convert some variables to factors and create new derived variables
data$Avoid_yn <- as.factor(ifelse(data$Resp_cat == "Avoid", 1, 0))
data$fem_repro <- as.factor(data$fem_repro)
data$Orientation_yn <- as.factor(data$Orientation_yn)

# Create a new variable for season (based on Pb_ID)
data <- data %>% 
  mutate(season = if_else(Pb_ID <= 18, "non-mating", "mating"))

# Standardize selected variables to have mean = 0, sd = 1
data_std <- data %>% 
  mutate(Familiarity = standardize(Familiarity),
         SRI = standardize(SRI),
         mal_strength = standardize(mal_strength),
         fem_strength = standardize(fem_strength),
         HR_overlap = standardize(HR_overlap),
         consort_dyad = standardize(consort_dyad),
         consort_male = standardize(consort_male),
         max_consort_length = standardize(max_consort_length),
         male_cumul_strength = standardize(male_cumul_strength),
         Resp_time = standardize(Resp_time))

#### Statistical Analyses ####

##### Model 1: Avoidance Distance Model #####

# Fit a linear mixed-effects model to predict avoidance distance
m_avoid <- lmer(Dist_avoid ~ Familiarity + male_cumul_strength + consort_dyad + SRI + HR_overlap + fem_strength + mal_strength + consort_male * fem_repro + (1 | Pb_signal_ID),
                data = data_std)

# Summarize model and check multicollinearity
summary(m_avoid)
vif(m_avoid)

# Model selection using dredge function
options(na.action = na.fail)
dr_m_avoid <- dredge(m_avoid)
options(na.action = na.omit)
top_m_avoid <- subset(dr_m_avoid, delta < 2)

# Model averaging for the top models
avg_m_avoid <- model.avg(top_m_avoid)
summary(avg_m_avoid)
confint(avg_m_avoid)

# Check for singularity issues in the model
isSingular(m_avoid)

# Check model assumptions: normality, heteroscedasticity, outliers
check_normality(m_avoid)
check_heteroscedasticity(m_avoid)
check_outliers(m_avoid)

# Visualize residuals vs fitted values
ggplot(data.frame(fitted = fitted(m_avoid), resid = residuals(m_avoid)), aes(fitted, resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed")

# QQ plot of residuals
qqPlot(residuals(m_avoid))

##### Model 2: Approach Distance Model #####

# Fit a generalized linear mixed model (GLMM) with a Tweedie distribution to predict approach distance
m_approach <- glmmTMB(Dist_approach ~ Familiarity + male_cumul_strength + consort_dyad + SRI + HR_overlap + fem_strength + mal_strength + consort_male * fem_repro + (1 | Pb_target) + (1 | Pb_signal_ID),
                      family = tweedie(link = "log"), 
                      data = data)

# Alternatively, using a linear mixed-effects model (if preferred)
m_approach <- lmer(Dist_approach ~ Familiarity + male_cumul_strength + consort_dyad + SRI + HR_overlap + fem_strength + mal_strength + consort_male * fem_repro + (1 | Pb_target) + (1 | Pb_signal_ID),
                   data = data_std)

# Summarize model and check multicollinearity
summary(m_approach)
vif(m_approach)

# Model selection using dredge function
options(na.action = na.fail)
dr_m_approach <- dredge(m_approach)
options(na.action = na.omit)

# Filter out models with non-converged results (AICc not available)
good_rows <- which(!is.na(dr_m_approach$AICc))
clean_dredge <- dr_m_approach[good_rows, ]

# Get top models from the clean subset
top_m_approach <- subset(clean_dredge, delta < 2)

# Model averaging for the top models
summary(model.avg(top_m_approach))
confint(model.avg(top_m_approach))

# Check for singularity issues in the model
isSingular(m_approach)

# Simulate residuals for validation
sim_res <- simulateResiduals(m_approach)
plot(sim_res)
testUniformity(sim_res)

##### Model 3: Duration of Response Model #####

# Fit a linear mixed-effects model for response duration
m_dur <- lmer(Resp_dur ~ Familiarity + male_cumul_strength + consort_male + Avoid_yn + consort_dyad + SRI + HR_overlap + fem_strength + mal_strength + consort_male * fem_repro + (1 | Pb_target) + (1 | Pb_signal_ID),
              data = data_std)

# Summarize model and check multicollinearity
summary(m_dur)
vif(m_dur)

# Model selection using dredge function
options(na.action = na.fail)
dr_m_dur <- dredge(m_dur)
options(na.action = na.omit)
top_m_dur <- subset(dr_m_dur, delta < 2)

# Model averaging for the top models
avg_m_dur <- model.avg(top_m_dur)

# Best model for response duration (simplified)
m_dur_best <- lmer(Resp_dur ~ Avoid_yn + fem_repro + (1 | Pb_target) + (1 | Pb_signal_ID),
                   data = data_std)
summary(m_dur_best)
confint(m_dur_best)

# Check for singularity issues in the model
isSingular(m_dur)

# Check model assumptions
check_normality(m_dur)
check_heteroscedasticity(m_dur)
check_outliers(m_dur)

# Visualize residuals vs fitted values
ggplot(data.frame(fitted = fitted(m_dur), resid = residuals(m_dur)), aes(fitted, resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed")

# QQ plot of residuals
qqPlot(residuals(m_dur))

##### Model 4: Orientation Response Model #####

# Fit a GLMM for orientation response
m_orient <- glmer(Orientation_yn ~ Familiarity + male_cumul_strength + consort_male + Avoid_yn + consort_dyad + SRI + HR_overlap + fem_strength + mal_strength + consort_male * fem_repro + (1 | Pb_target) + (1 | Pb_signal_ID),
                  data = data_std, family = binomial())

# Summarize model
summary(m_orient)

# Model selection using dredge function
options(na.action = na.fail)
dr_m_orient <- dredge(m_orient)
options(na.action = na.omit)
top_m_orient <- subset(dr_m_orient, delta < 2)

# Model averaging for the top models
summary(model.avg(top_m_orient))
confint(model.avg(top_m_orient))

# Check for singularity issues in the model
isSingular(m_orient)

# Simulate residuals for model validation
plot(simulateResiduals(m_orient))

# Check for overdispersion
check_overdispersion(m_orient)
