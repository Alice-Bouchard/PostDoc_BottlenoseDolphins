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

# Model selection using dredge function
options(na.action = na.fail)
dr_m_avoid <- dredge(m_avoid)
options(na.action = na.omit)
top_m_avoid <- subset(dr_m_avoid, delta < 2)

# Model averaging for the top models
avg_m_avoid <- model.avg(top_m_avoid)
summary(avg_m_avoid)
confint(avg_m_avoid)

##### Model 2: Approach Distance Model #####

# Fit a generalized linear mixed model (GLMM) with a Tweedie distribution to predict approach distance
m_approach <- glmmTMB(Dist_approach ~ Familiarity + male_cumul_strength + consort_dyad + SRI + HR_overlap + fem_strength + mal_strength + consort_male * fem_repro + (1 | Pb_target) + (1 | Pb_signal_ID),
                      family = tweedie(link = "log"), 
                      data = data)

# Summarize model and check multicollinearity
summary(m_approach)
vif(m_approach)

# Check for singularity issues in the model
isSingular(m_approach)

# Simulate residuals for validation
sim_res <- simulateResiduals(m_approach)
plot(sim_res)
testUniformity(sim_res)

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

##### Model 3: Duration of Response Model #####

# Fit a linear mixed-effects model for response duration
m_dur <- lmer(Resp_dur ~ Familiarity + male_cumul_strength + consort_male + Avoid_yn + consort_dyad + SRI + HR_overlap + fem_strength + mal_strength + consort_male * fem_repro + (1 | Pb_target) + (1 | Pb_signal_ID),
              data = data_std)

# Summarize model and check multicollinearity
summary(m_dur)
vif(m_dur)

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

# Model selection using dredge function
options(na.action = na.fail)
dr_m_dur <- dredge(m_dur)
options(na.action = na.omit)
top_m_dur <- subset(dr_m_dur, delta < 2)

# Summarize top model for response duration (only one model in top models so no model averaging)
m_dur_best <- lmer(Resp_dur ~ Avoid_yn + fem_repro + (1 | Pb_target) + (1 | Pb_signal_ID),
                   data = data_std)
summary(m_dur_best)
confint(m_dur_best)


##### Model 4: Orientation Response Model #####

# Fit a GLMM for orientation response
m_orient <- glmer(Orientation_yn ~ Familiarity + male_cumul_strength + consort_male + Avoid_yn + consort_dyad + SRI + HR_overlap + fem_strength + mal_strength + consort_male * fem_repro + (1 | Pb_target) + (1 | Pb_signal_ID),
                  data = data_std, family = binomial())

# Summarize model and check multicollinearity
summary(m_orient)
vif(m_orient)

# Check for singularity issues in the model
isSingular(m_orient)

# Simulate residuals for model validation
plot(simulateResiduals(m_orient))

# Check for overdispersion
check_overdispersion(m_orient)

# Model selection using dredge function
options(na.action = na.fail)
dr_m_orient <- dredge(m_orient)
options(na.action = na.omit)
top_m_orient <- subset(dr_m_orient, delta < 2)

# Model averaging for the top models
summary(model.avg(top_m_orient))
confint(model.avg(top_m_orient))

#### graphs ####

# Graph 1: Avoidance distance vs. Male consortship rate, annotated with p-value
ggplot(data, aes(x = consort_male, y = Dist_avoid, col = as.factor(fem_repro))) + 
  geom_point() + 
  geom_smooth(method = 'lm', aes(fill = as.factor(fem_repro)), alpha = 0.1) +
  labs(
    y = "Avoidance distance (m)", 
    x = "Male consortship rate", 
    col = "Female reproductive state"
  ) +
  scale_color_manual(values = c("#00B89F", "#D6743E"), labels = c("Non available", "Available")) +
  scale_fill_manual(values = c("#00B89F", "#D6743E")) +
  annotate("text", x = 0.47, y = 150, label = "p = 0.002", col = "#D6743E") +
  guides(fill = "none") +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, size = 1),
    legend.position = c(0.3, 0.85)
  )

# Graph 2: Response duration vs. Female reproductive state (Boxplot)
ggplot(data, aes(x = fem_repro, y = Resp_dur, col = fem_repro)) +
  geom_boxplot() + 
  geom_jitter(alpha = 0.3) +
  labs(
    y = "Response duration (s)", 
    x = "Female reproductive state"
  ) +
  scale_color_manual(values = c("#00B89F", "#D6743E")) +
  scale_x_discrete(labels = c("Not available", "Available")) +
  annotate("text", x = 2, y = 98, label = "p = 0.005", col = "#D6743E") +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, size = 1),
    legend.position = "none"
  )

# Graph 3: Response duration vs. Avoidance response (Boxplot)
ggplot(data, aes(x = Avoid_yn, y = Resp_dur, col = Avoid_yn)) +
  geom_boxplot() + 
  geom_jitter(alpha = 0.3) +
  labs(
    y = "Response duration (s)", 
    x = "Avoidance response"
  ) +
  scale_color_manual(values = c("#124E78", "#D03970")) +
  scale_x_discrete(labels = c("No", "Yes")) +
  annotate("text", x = 2, y = 98, label = "p = 0.008", col = "#D03970") +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, size = 1),
    legend.position = "none"
  )

# Graph 4a: Male ID vs. Approach distance (Boxplot and Jitter)
plot_male_approach <- ggplot(data, aes(x = Pb_signal_ID, y = Dist_approach, col = Pb_signal_ID)) + 
  geom_boxplot() + 
  geom_jitter(alpha = 0.4) +
  scale_color_manual(values = paletteer_c("ggthemes::Orange-Blue Diverging", 11)) +
  xlab("Male ID") + 
  ylab("Approach distance (m)") +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, size = 1),
    legend.position = "none"
  )

# Graph 4b: Male ID vs. Avoidance distance (Boxplot and Jitter)
plot_male_avoid <- ggplot(data, aes(x = Pb_signal_ID, y = Dist_avoid, col = Pb_signal_ID)) + 
  geom_boxplot() + 
  geom_jitter(alpha = 0.4) +
  xlab(NULL) + 
  ylab("Avoidance distance (m)") +
  scale_color_manual(values = paletteer_c("ggthemes::Orange-Blue Diverging", 11)) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, size = 1),
    legend.position = "none"
  )

# Combine the Male ID graphs into a single figure
ggarrange(
  plot_male_avoid, 
  plot_male_approach, 
  nrow = 2, 
  labels = c("(a)", "(b)"),
  label.x = -0.01,
  font.label = list(size = 12, color = "black", face = "bold", family = NULL)
)

