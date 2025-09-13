library(nlme)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(lme4)
library(effectsize)
library(psych)
library(tidyr)
library(dplyr)


importDat <- read.csv('/Users/yanbinniu/local_project/social_flanker_eeg/script/analyses_erp/data/beh_all_0224.csv', header = TRUE)

##################### 1) ACC ##################### 
# 1. long formatting
# Organize into long format for ACC: subject, acc, congruency, observation
# "subject" "acc_congruent_nonsocial" "acc_incongruent_nonsocial" "acc_congruent_social" "acc_incongruent_social"
# Step 1: Subset the relevant columns
subsetDat <- importDat %>%
  select(subject, acc_congruent_nonsocial, acc_incongruent_nonsocial, acc_congruent_social, acc_incongruent_social)

# Step 2: Convert to long format
longDat_acc <- subsetDat %>%
  pivot_longer(cols = starts_with("acc_"),  # Select columns that start with "acc_"
               names_to = c("congruency", "observation"),  # Create new columns from the parts of names
               names_pattern = "acc_(.*)_(.*)",  # Regex pattern to extract congruency and observation from column names
               values_to = "acc")  # Name for the column holding the values

# 2. contrast and standardize
longDat_acc$subject <- as.factor(longDat_acc$subject)
longDat_acc$congruency <- as.factor(longDat_acc$congruency)
longDat_acc$observation <- as.factor(longDat_acc$observation)
contrasts(longDat_acc$congruency) <- rev(contr.sum(2))
contrasts(longDat_acc$observation) <- rev(contr.sum(2))
longDat_acc$acc_z <- scale(longDat_acc$acc , center = TRUE, scale = TRUE)

# 3. modeling
fit_model <- lme(acc_z ~ congruency*observation, 
                 random = ~ 1 | subject, 
                 data=longDat_acc, 
                 method="ML")
summary(fit_model)
intervals(fit_model)

# 4. follow up
# 4.1 follow up main effect of congruency
library(emmeans)
fit_mean <- lme(acc ~ congruency*observation, 
                random = ~ 1 | subject, 
                data=longDat_acc, 
                method="ML")
# Compute the estimated marginal means (mean accuracy) for the observation variable
emms <- emmeans(fit_mean, ~ congruency, type = "response")
# Display the estimated marginal means
print(emms)
# congruency  emmean      SE df lower.CL upper.CL
# congruent    0.949 0.00918 53    0.930    0.967
# incongruent  0.795 0.00918 53    0.776    0.813
# Format the emmeans output with 4 digits
emms_df <- as.data.frame(summary(emms))
format(emms_df, digits = 4, nsmall = 4)

##################### 2) RT ##################### 
importDat <- read.csv('/Users/yanbinniu/local_project/social_flanker_eeg/script/analyses_erp/data/beh_all_0224.csv', header = TRUE)

# 1. long formatting
# Step 1: Subset the relevant columns
subsetDat <- importDat %>%
  select(subject, 
         rt_congruent_correct_nonsocial, rt_congruent_error_nonsocial, 
         rt_incongruent_correct_nonsocial, rt_incongruent_error_nonsocial, 
         rt_congruent_correct_social, rt_congruent_error_social, 
         rt_incongruent_correct_social, rt_incongruent_error_social)

# Step 2: Convert to long format
longDat_rt <- subsetDat %>%
  pivot_longer(cols = starts_with("rt_"),  # Select columns that start with "rt_"
               names_to = c("congruency", "accuracy", "observation"),  # Create new columns from parts of the names
               names_pattern = "rt_(.*)_(.*)_(.*)",  # Regex pattern to extract congruency, accuracy, and observation
               values_to = "rt")  # Name for the column holding the values

# 2. contrast and standardize
longDat_rt$subject <- as.factor(longDat_rt$subject)
longDat_rt$congruency <- as.factor(longDat_rt$congruency)
longDat_rt$observation <- as.factor(longDat_rt$observation)
longDat_rt$accuracy <- as.factor(longDat_rt$accuracy)
contrasts(longDat_rt$congruency) <- rev(contr.sum(2))
contrasts(longDat_rt$observation) <- rev(contr.sum(2))
contrasts(longDat_rt$accuracy) <- rev(contr.sum(2))
longDat_rt$log_rt <- log(longDat_rt$rt*1000 + 1)
longDat_rt$log_rt_z <- scale(longDat_rt$log_rt , center = TRUE, scale = TRUE)

# 3. modeling
fit_model <- lme(log_rt_z ~ congruency*observation*accuracy, 
                 random = ~ 1 | subject, 
                 data=longDat_rt, 
                 method="ML", 
                 na.action = na.omit)
summary(fit_model)
intervals(fit_model)

# 4. follow up
# 4.1 follow up main effect of congruency
library(emmeans)
fit_mean <- lme(log_rt ~ congruency*observation*accuracy, 
                random = ~ 1 | subject, 
                data=longDat_rt, 
                method="ML", 
                na.action = na.omit)
# Compute the estimated marginal means (mean accuracy) for the observation variable
emms <- emmeans(fit_mean, ~ congruency, type = "response")
# Display the estimated marginal means
print(emms)
# congruency  emmean     SE df lower.CL upper.CL
# congruent     5.92 0.0246 53     5.87     5.97
# incongruent   6.03 0.0244 53     5.98     6.08
# > exp(5.92)
# [1] 372.4117
# > exp(6.03)
# [1] 415.715

# 4.2 follow up main effect of accuracy
# Compute the estimated marginal means (mean accuracy) for the observation variable
emms <- emmeans(fit_mean, ~ accuracy, type = "response")
# Display the estimated marginal means
print(emms)
# accuracy emmean     SE df lower.CL upper.CL
# correct    6.06 0.0244 53     6.01     6.11
# error      5.89 0.0246 53     5.84     5.94
# > exp(6.06)
# [1] 428.3754
# > exp(5.89)
# [1] 361.4053

# 4.3 follow up congruency*accuracy interaction
# Compute the estimated marginal means (mean accuracy) for the observation variable
emms <- emmeans(fit_mean, ~ congruency*accuracy, type = "response")
# Display the estimated marginal means
print(emms)
# congruency  accuracy emmean     SE df lower.CL upper.CL
# congruent   correct    5.97 0.0265 53     5.92     6.02
# incongruent correct    6.14 0.0265 53     6.09     6.20
# congruent   error      5.86 0.0270 53     5.81     5.92
# incongruent error      5.92 0.0265 53     5.86     5.97
# [1] 464.0536
# > exp(5.97)
# [1] 391.5057
# > exp(6.14)
# [1] 464.0536
# > exp(5.86)
# [1] 350.7241
# > exp(5.92)
# [1] 372.4117
# Pairwise comparisons of slopes
slopes_comparisons <- pairs(emms)
summary(slopes_comparisons)
# contrast                                estimate     SE  df t.ratio p.value
# congruent correct - incongruent correct  -0.1761 0.0207 361  -8.512  <.0001
# congruent correct - congruent error       0.1053 0.0213 361   4.951  <.0001
# congruent correct - incongruent error     0.0509 0.0207 361   2.461  0.0679
# incongruent correct - congruent error     0.2814 0.0213 361  13.231  <.0001
# incongruent correct - incongruent error   0.2270 0.0207 361  10.974  <.0001
# congruent error - incongruent error      -0.0544 0.0213 361  -2.557  0.0532

