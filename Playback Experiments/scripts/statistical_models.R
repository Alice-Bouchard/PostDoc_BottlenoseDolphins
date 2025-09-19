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
library(datawizard) # Data transformations (standardization)
library(DHARMa)     # Residual diagnostics for mixed models
library(glmmTMB)    # Generalized linear mixed-effects models
library(ggeffects)  # Marginal effects and predicted values with ggplot2 support

# Load the dataset
setwd("~/Perso/PostDoc")
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
m_avoid <- lmer(Dist_avoid ~ Familiarity + male_cumul_strength + consort_dyad + SRI* fem_repro + HR_overlap + fem_strength + mal_strength + consort_male * fem_repro + (1 | Pb_target) + (1 | Pb_signal_ID),
                data = data_std)

# Summarize model and check multicollinearity
summary(m_avoid)
vif(m_avoid)

# Check for singularity issues in the model
isSingular(m_avoid)
# Check model
qqPlot(residuals(m_avoid))
plot(m_avoid)


# Model selection using dredge function
options(na.action = na.fail)
dr_m_avoid <- dredge(m_avoid)
options(na.action = na.omit)
top_m_avoid <- subset(dr_m_avoid, delta < 2)

# Model averaging for the top models
summary(model.avg(top_m_avoid))
confint(model.avg(top_m_avoid))

##### Model 2: Approach Distance Model #####

# Fit a generalized linear mixed model (GLMM) with a Tweedie distribution to predict approach distance
m_approach <- glmmTMB(Dist_approach ~ Familiarity + male_cumul_strength + consort_dyad + SRI * fem_repro  + HR_overlap + fem_strength + mal_strength + consort_male * fem_repro + (1 | Pb_target) + (1 | Pb_signal_ID),
                      family = tweedie(link = "log"), 
                      data = data_std)

# Summarize model 
summary(m_approach)

# Simulate residuals for validation
sim_res <- simulateResiduals(m_approach)
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
m_dur <- lmer(Resp_dur ~ Familiarity + male_cumul_strength + consort_male + Avoid_yn + consort_dyad + SRI* fem_repro + HR_overlap + fem_strength + mal_strength + consort_male * fem_repro  + (1 | Pb_target) + (1 | Pb_signal_ID),
              data = data_std)

# Summarize model and check multicollinearity
summary(m_dur)
vif(m_dur)

# Check for singularity issues in the model
isSingular(m_dur)
# Check model
qqPlot(residuals(m_dur))
plot(m_dur)

# Model selection using dredge function
options(na.action = na.fail)
dr_m_dur <- dredge(m_dur)
options(na.action = na.omit)
top_m_dur <- subset(dr_m_dur, delta < 2)

# Model averaging for the top models
summary(model.avg(top_m_dur))
confint(model.avg(top_m_dur))


##### Model 4: Orientation Response Model #####

# Fit a GLMM for orientation response
m_orient <- glmer(Orientation_yn ~ Familiarity + male_cumul_strength + consort_male + Avoid_yn + consort_dyad + SRI * fem_repro + HR_overlap + fem_strength + mal_strength + consort_male * fem_repro + (1 | Pb_target) + (1 | Pb_signal_ID),
                  data = data_std, family = binomial())

# Summarize model and check multicollinearity
summary(m_orient)
vif(m_orient)

# Check for singularity issues in the model
isSingular(m_orient)
# Check model
plot(simulateResiduals(m_orient))# Simulate residuals for model validation
check_overdispersion(m_orient)# Check for overdispersion

# Model selection using dredge function
options(na.action = na.fail)
dr_m_orient <- dredge(m_orient)
options(na.action = na.omit)
top_m_orient <- subset(dr_m_orient, delta < 2)

# Model averaging for the top models
summary(model.avg(top_m_orient))
confint(model.avg(top_m_orient))

#### graphs ####

# Graph 1: Avoidance distance vs. Male consortship rate

# Create a newdata grid for predictions
newdata1 <- expand.grid(
  consort_male = seq(min(data_std$consort_male), max(data_std$consort_male), length.out = 100),
  fem_repro = factor(c(0, 1), levels = levels(data_std$fem_repro)),
  Familiarity = mean(data_std$Familiarity, na.rm = TRUE),
  male_cumul_strength = mean(data_std$male_cumul_strength, na.rm = TRUE),
  consort_dyad = mean(data_std$consort_dyad, na.rm = TRUE),
  HR_overlap = mean(data_std$HR_overlap, na.rm = TRUE),
  fem_strength = mean(data_std$fem_strength, na.rm = TRUE),
  mal_strength = mean(data_std$mal_strength, na.rm = TRUE),
  SRI = mean(data_std$SRI, na.rm = TRUE)
)

# Get predicted values from the model (exclude random effects)
preds1 <- as.data.frame(predict(m_avoid, newdata = newdata1, re.form = NA, se.fit = TRUE))

# Combine predictions with newdata
preds1 <- cbind(newdata1, preds1)

# Extract mean and sd from original variable
mean_c <- mean(data$consort_male, na.rm = TRUE)
sd_c   <- sd(data$consort_male, na.rm = TRUE)

# Back-transform for plotting
preds1$consort_male_original <- preds1$consort_male * sd_c + mean_c

# Plot raw points + predicted slopes
plot1 <- ggplot() +
  geom_point(data = data, 
             aes(x = consort_male, y = Dist_avoid, colour = factor(fem_repro)), alpha = 0.6) +
  geom_line(data = preds1, 
            aes(x = consort_male_original, y = fit, colour = factor(fem_repro)), size = 0.8) +
  geom_ribbon(data = preds1, 
              aes(x = consort_male_original, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit, fill = factor(fem_repro)),
              alpha = 0.1, colour = NA) +
  annotate("text", x = 0.45, y = 180, label = "p < 0.001", col = "#D6743E") +
  labs(x = "Male consortship rate", y = "Avoidance distance (m)", colour = "Reproductive state", fill = "Reproductive state") +
  scale_color_manual(values = c("#00B89F", "#D6743E"), labels = c("Non available", "Available")) +
  scale_fill_manual(values = c("#00B89F", "#D6743E")) +
  guides(fill = "none") +
  labs(tag = "(a)") +
  scale_y_continuous(breaks = function(x) seq(0, max(x), by = 100)) +
  theme_bw()+ 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, size = 1),
    legend.position = c(0.2, 0.75)
  )

# Graph 2: Avoidance distance vs. SRI

# Create a newdata grid for predictions
newdata2 <- expand.grid(
  SRI = seq(min(data_std$SRI), max(data_std$SRI), length.out = 100),
  fem_repro = factor(c(0, 1), levels = levels(data_std$fem_repro)),  
  Familiarity = mean(data_std$Familiarity, na.rm = TRUE),
  male_cumul_strength = mean(data_std$male_cumul_strength, na.rm = TRUE),
  consort_dyad = mean(data_std$consort_dyad, na.rm = TRUE),
  HR_overlap = mean(data_std$HR_overlap, na.rm = TRUE),
  fem_strength = mean(data_std$fem_strength, na.rm = TRUE),
  mal_strength = mean(data_std$mal_strength, na.rm = TRUE),
  consort_male = mean(data_std$consort_male, na.rm = TRUE)
)

# Get predicted values from the model (exclude random effects)
preds2 <- as.data.frame(predict(m_avoid, newdata = newdata2, re.form = NA, se.fit = TRUE))

# Combine predictions with newdata
preds2 <- cbind(newdata2, preds2)

# Extract mean and sd from original variable
mean_s <- mean(data$SRI, na.rm = TRUE)
sd_s   <- sd(data$SRI, na.rm = TRUE)

# Back-transform for plotting
preds2$SRI_original <- preds2$SRI * sd_s + mean_s


# Plot raw points + predicted slopes
plot2 <- ggplot() +
  geom_point(data = data, 
             aes(x = SRI, y = Dist_avoid, colour = factor(fem_repro)), alpha = 0.6) +
  geom_line(data = preds2, 
            aes(x = SRI_original, y = fit, colour = factor(fem_repro)), size = 0.8) +
  geom_ribbon(data = preds2, 
              aes(x = SRI_original, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit, fill = factor(fem_repro)),
              alpha = 0.1, colour = NA) +
  annotate("text", x = 0.07, y = 265, label = "p = 0.004", col = "#D6743E") +
  labs(x = "SRI", y = "Avoidance distance (m)", colour = "Reproductive state", fill = "Reproductive state") +
  scale_color_manual(values = c("#00B89F", "#D6743E"), labels = c("Non available", "Available")) +
  scale_fill_manual(values = c("#00B89F", "#D6743E")) +
  guides(fill = "none") +
  theme_bw()+ 
  labs(tag = "(b)") +
  scale_x_continuous(
    breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10)  # only these values will have labels
  ) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, size = 1),
    legend.position = c(0.2, 0.75)
  )

# Graph 3: Response duration vs. Male consortship rate

# Create a newdata grid for predictions
newdata3 <- expand.grid(
  consort_male = seq(min(data_std$consort_male), max(data_std$consort_male), length.out = 100),
  fem_repro = factor(c(0, 1), levels = levels(data_std$fem_repro)),
  Avoid_yn = factor("0", levels = levels(data_std$Avoid_yn)),
  Familiarity = mean(data_std$Familiarity, na.rm = TRUE),
  male_cumul_strength = mean(data_std$male_cumul_strength, na.rm = TRUE),
  consort_dyad = mean(data_std$consort_dyad, na.rm = TRUE),
  HR_overlap = mean(data_std$HR_overlap, na.rm = TRUE),
  fem_strength = mean(data_std$fem_strength, na.rm = TRUE),
  mal_strength = mean(data_std$mal_strength, na.rm = TRUE),
  SRI = mean(data_std$SRI, na.rm = TRUE)
)

# Get predicted values from the model (exclude random effects)
preds3 <- as.data.frame(predict(m_dur, newdata = newdata3, re.form = NA, se.fit = TRUE))

# Combine predictions with newdata
preds3 <- cbind(newdata3, preds3) 

# Back-transform for plotting
preds3$consort_male_original <- preds3$consort_male * sd_c + mean_c

# Plot raw points + predicted slopes
plot3 <- ggplot() +
  geom_point(data = data, 
             aes(x = consort_male, y = Resp_dur, colour = factor(fem_repro)), alpha = 0.6) +
  geom_line(data = preds3, 
            aes(x = consort_male_original, y = fit, colour = factor(fem_repro), group = fem_repro), size = 0.8) +
  geom_ribbon(data = preds3, 
              aes(x = consort_male_original, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit, fill = factor(fem_repro), group = fem_repro),
              alpha = 0.1, colour = NA) +
  annotate("text", x = 0.45, y = 110, label = "p = 0.01", col = "#D6743E") +
  labs(x = "Male consortship rate", y = "Response duration (s)", colour = "Reproductive state", fill = "Reproductive state") +
  scale_color_manual(values = c("#00B89F", "#D6743E"), labels = c("Non available", "Available")) +
  scale_fill_manual(values = c("#00B89F", "#D6743E")) +
  guides(fill = "none") +
  theme_bw()+ 
  labs(tag = "(c)") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, size = 1),
    legend.position = c(0.2, 0.75)
  )

# Combine the graphs into a single figure
fig1 <- ggarrange(plot1,plot2,plot3,nrow = 3)

# Save the figure in png format
ggsave("fig1.png", plot = fig1, width = 195, height = 240, units = "mm", dpi = 300)

# Graph on male ID

# Reorder Pb_signal_ID factor levels by median avoidance distance
data <- data %>%
  mutate(Pb_signal_ID = reorder(Pb_signal_ID, Dist_avoid, FUN = median))


# Graph 4a: Male ID vs. Avoidance distance (Boxplot and Jitter)
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

# Graph 4b: Male ID vs. Approach distance (Boxplot and Jitter)
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

# Combine the Male ID graphs into a single figure
fig2 <- ggarrange(
  plot_male_avoid, 
  plot_male_approach, 
  nrow = 2, 
  labels = c("(a)", "(b)"),
  label.x = -0.01,
  font.label = list(size = 12, color = "black", face = "bold", family = NULL)
)

# Save the figure in png format
ggsave("fig2.png", plot = fig2, width = 195, height = 150, units = "mm", dpi = 300)

