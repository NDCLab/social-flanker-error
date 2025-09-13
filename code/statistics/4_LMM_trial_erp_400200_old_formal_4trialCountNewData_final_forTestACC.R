library(nlme)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(lme4)
library(effectsize)
library(psych)
library(tidyr)
library(dplyr)


importDat <- read.csv('/Users/yanbinniu/local_project/social_flanker_eeg/script/analyses_erp/data/df_fne_erp_merged_oldChn_400And200_with_trialNO.csv', header = TRUE)
dat4Analyses <- importDat[which(importDat$congruency =='i' & importDat$validTrial == 1), ]

# Factorize the categorical variables
dat4Analyses$subject <- as.factor(dat4Analyses$subject)
dat4Analyses$observation <- as.factor(dat4Analyses$observation)
dat4Analyses$accuracy <- as.factor(dat4Analyses$accuracy)

# Assign the sum contrast codes
contrasts(dat4Analyses$observation) <- rev(contr.sum(2))
contrasts(dat4Analyses$accuracy) <- rev(contr.sum(2))

# 1. check if accuracy*observation matters for trial counts
table_incongruent <- table(dat4Analyses$subject,dat4Analyses$accuracy, dat4Analyses$observation)
write.table(table_incongruent, file = "incongruent_trial_count_new.csv", sep = ",", quote = FALSE, row.names = F)
incongruent_count <- read.csv('incongruent_trial_count_new.csv', header = TRUE)
# assigning new names to the columns of the data frame
colnames(incongruent_count) <- c('subj','acc','observation','Freq')
incongruent_count$acc <- as.factor(incongruent_count$acc)
incongruent_count$observation <- as.factor(incongruent_count$observation)
# fit
fit_count <- lme(Freq ~ observation*acc, random = ~ 1 | subj, data=incongruent_count)
anova(fit_count)
summary(fit_count)
plot_model(fit_count, type = "int", terms = c("acc", "observation"))
# subject to remove -- less than 6
# subj_to_remove <- unique(incongruent_count[which((incongruent_count$Freq<6)&(incongruent_count$Freq>0)),]$subj)
cond_to_remove <- incongruent_count[which((incongruent_count$Freq<6)&(incongruent_count$Freq>0)),]
cond_to_remove
#       subj acc observation Freq
# 4   160018   0          ns    4
# 9   160023   0          ns    2
# 38  160068   0          ns    3
# 45  160076   0          ns    2
# 71  160049   1          ns    4
# 97  160021   0           s    5
# 103 160029   0           s    4
# 106 160033   0           s    5
# 124 160063   0           s    2

# 3. remove subj with less 6 trials: 49 - 14 = 35 subjects
dat4Analyses_reducedby_8 <- dat4Analyses[!(((dat4Analyses$subject=='160018') & (dat4Analyses$observation=='ns'))
                                           | ((dat4Analyses$subject=='160023') & (dat4Analyses$observation=='ns'))
                                           | ((dat4Analyses$subject=='160068') & (dat4Analyses$observation=='ns'))
                                           | ((dat4Analyses$subject=='160076') & (dat4Analyses$observation=='ns'))
                                           | ((dat4Analyses$subject=='160049') & (dat4Analyses$observation=='ns'))
                                           | ((dat4Analyses$subject=='160021') & (dat4Analyses$observation=='s'))
                                           | ((dat4Analyses$subject=='160029') & (dat4Analyses$observation=='s'))
                                           | ((dat4Analyses$subject=='160033') & (dat4Analyses$observation=='s'))
                                           | ((dat4Analyses$subject=='160063') & (dat4Analyses$observation=='s'))),]

# contrast and standardize
contrasts(dat4Analyses_reducedby_8$observation) <- rev(contr.sum(2))
contrasts(dat4Analyses_reducedby_8$accuracy) <- rev(contr.sum(2))
dat4Analyses_reducedby_8$erp_ern <- scale(dat4Analyses_reducedby_8$erp_ern , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$age_yr <- scale(dat4Analyses_reducedby_8$age_yr , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$bfne_sum <- scale(dat4Analyses_reducedby_8$bfne_sum , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$scaared_social_phobia <- scale(dat4Analyses_reducedby_8$scaared_social_phobia , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$stai_sum <- scale(dat4Analyses_reducedby_8$stai_sum , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$stai_sum2 <- scale(dat4Analyses_reducedby_8$stai_sum2 , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$stai_mean <- scale(dat4Analyses_reducedby_8$stai_mean , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$num_all <- scale(dat4Analyses_reducedby_8$num_all , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$num_block <- scale(dat4Analyses_reducedby_8$num_block , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$num_inblock <- scale(dat4Analyses_reducedby_8$num_inblock , center = TRUE, scale = TRUE)


#########################################################################
# 1. run main model

#########################################################################
# 1) main model
fit_dat4Analyses_reducedby_8 <- glmer(accuracy ~ observation*num_block*num_inblock + (1 | subject), 
                                      data = dat4Analyses_reducedby_8, 
                                      family = binomial(link = "logit"),
                                      na.action = na.omit)
# (1 | subject/num_block)
summary(fit_dat4Analyses_reducedby_8)
# Get Wald-based confidence intervals
conf_intervals_wald <- confint(fit_dat4Analyses_reducedby_8, method = "Wald")
print(conf_intervals_wald)

# 2) follow up  main effect of observation
library(emmeans)
# Compute the estimated marginal means (mean accuracy) for the observation variable
emms <- emmeans(fit_dat4Analyses_reducedby_8, ~ observation, type = "response")
# Display the estimated marginal means
print(emms)

