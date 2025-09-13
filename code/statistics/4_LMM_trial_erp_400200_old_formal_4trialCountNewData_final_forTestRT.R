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
# write.table(table_incongruent, file = "incongruent_trial_count_new.csv", sep = ",", quote = FALSE, row.names = F)
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
dat4Analyses$observation <- as.factor(dat4Analyses$observation)
dat4Analyses$accuracy <- as.factor(dat4Analyses$accuracy)
contrasts(dat4Analyses_reducedby_8$observation) <- rev(contr.sum(2))
contrasts(dat4Analyses_reducedby_8$accuracy) <- rev(contr.sum(2))
# use log(rt + 1) to avoid undefined log(0), change to ms first
dat4Analyses_reducedby_8$log_rt <- log(dat4Analyses_reducedby_8$rt*1000 + 1)
dat4Analyses_reducedby_8$log_rt_z <- scale(dat4Analyses_reducedby_8$log_rt , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$num_all_z <- scale(dat4Analyses_reducedby_8$num_all , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$num_block_z <- scale(dat4Analyses_reducedby_8$num_block , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$num_inblock_z <- scale(dat4Analyses_reducedby_8$num_inblock , center = TRUE, scale = TRUE)


#########################################################################
# 1. run main model

#########################################################################
# 1) main model
fit_dat4Analyses_reducedby_8 <- lme(log_rt_z ~ accuracy*observation*num_block_z*num_inblock_z, 
                                    random = ~ 1 | subject, 
                                    data=dat4Analyses_reducedby_8, 
                                    method="ML")
summary(fit_dat4Analyses_reducedby_8)
intervals(fit_dat4Analyses_reducedby_8)
reghelper::simple_slopes(fit_dat4Analyses_reducedby_8, 
                         levels = list(num_block = 'sstest',
                                       accuracy = c(0, 1)),
                         confint = T)

#####################
# 2) Follow Up Analysis
#####################
# 1. main effect of accuracy
library(emmeans)
fit_mean <- lme(log_rt ~ accuracy*observation*num_block*num_inblock, 
                                    random = ~ 1 | subject, 
                                    data=dat4Analyses_reducedby_8, 
                                    method="ML")
# Compute the estimated marginal means (mean accuracy) for the observation variable
emms <- emmeans(fit_mean, ~ accuracy, type = "response")
# Display the estimated marginal means
print(emms)
# accuracy emmean     SE df lower.CL upper.CL
# 0          5.90 0.0248 43     5.85     5.95
# 1          6.14 0.0244 43     6.09     6.19
# > exp(5.90)
# [1] 365.0375
# > exp(6.14)
# [1] 464.0536

#####################
# 2 -- slopes for num_block*accuracy
# Estimate slopes for num_block at accuracy
# setting 'at' for fixed value of num_block and calculating slopes for num_inblock
trends <- emtrends(fit_dat4Analyses_reducedby_8, ~ accuracy, 
                   "num_block_z", 
                   at = list(num_inblock_z = mean(dat4Analyses_reducedby_8$num_inblock_z)))
summary(trends)
tests <- test(trends, null = 0)
summary(tests)

# contrasts
# Set up the emmeans for the interaction of interest, specifically looking at slopes
mean_num_inblock_z <- mean(dat4Analyses_reducedby_8$num_inblock_z, na.rm = TRUE)
emm_slopes <- emtrends(fit_dat4Analyses_reducedby_8, 
                       specs = pairwise ~ accuracy, 
                       var = "num_block_z", 
                       at = list(num_inblock_z = mean_num_inblock_z))
summary(emm_slopes)

# Pairwise comparisons of slopes
slopes_comparisons <- pairs(emm_slopes)
summary(slopes_comparisons)
# Adjust for multiple testing if necessary
summary(slopes_comparisons, adjust = "bonferroni")

#####################
# 3 -- slopes for num_block*observation
# Estimate slopes for num_block at observation
# setting 'at' for fixed value of num_block and calculating slopes for num_inblock
trends <- emtrends(fit_dat4Analyses_reducedby_8, ~ observation, 
                   "num_block_z", 
                   at = list(num_inblock_z = mean(dat4Analyses_reducedby_8$num_inblock_z)))
summary(trends)
tests <- test(trends, null = 0)
summary(tests)

