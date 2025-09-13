library(nlme)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(lme4)
library(effectsize)
library(psych)
library(tidyr)
library(dplyr)

##################### old code for data prepartion #####################
# # switch to the desired workng directory
# setwd('/Volumes/T7_2T/mac_backup_02262024/Projects/social_flanker_eeg/script/analyses_erp')

# importDat <- read.csv('/Volumes/T7_2T/mac_backup_02262024/Projects/social_flanker_eeg/script/analyses_erp/data/df_fne_erp_merged_oldChn_400And200_with_trialNO_fd.csv', header = TRUE)
importDat <- read.csv('/Users/yanbinniu/local_project/social_flanker_eeg/script/analyses_erp/data/df_fne_erp_merged_oldChn_400And200_with_trialNO_fd.csv', header = TRUE)
dat4Analyses <- importDat[which(importDat$congruency =='i' & importDat$validTrial == 1), ]

# Factorize the categorical variables
dat4Analyses$subject <- as.factor(dat4Analyses$subject)
dat4Analyses$observation <- as.factor(dat4Analyses$observation)
dat4Analyses$accuracy <- as.factor(dat4Analyses$accuracy)

# Assign the sum contrast codes
contrasts(dat4Analyses$observation) <- rev(contr.sum(2))
contrasts(dat4Analyses$accuracy) <- rev(contr.sum(2))

# 1. check if accuracy*observation affects trial counts
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
dat4Analyses_reducedby_8$erp_ern_z <- scale(dat4Analyses_reducedby_8$erp_ern, center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$num_all_z <- scale(dat4Analyses_reducedby_8$num_all , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$num_block_z <- scale(dat4Analyses_reducedby_8$num_block , center = TRUE, scale = TRUE)
dat4Analyses_reducedby_8$num_inblock_z <- scale(dat4Analyses_reducedby_8$num_inblock , center = TRUE, scale = TRUE)

# Create the new column 'go_faster_contrast' based on 'feedback'
dat4Analyses_reducedby_8$go_faster_contrast <- ifelse(dat4Analyses_reducedby_8$feedback == 1, -1,
                                                      ifelse(dat4Analyses_reducedby_8$feedback == 2, 1, NA))
dat4Analyses_reducedby_8$go_faster_contrast <- factor(dat4Analyses_reducedby_8$go_faster_contrast, levels = c(-1, 1))
contrasts(dat4Analyses_reducedby_8$go_faster_contrast) <- rev(contr.sum(2))


#########################################################################
# 1. run main model

#########################################################################
# 1) main model
fit_dat4Analyses_reducedby_8 <- lme(erp_ern_z ~ accuracy*observation*num_block_z*num_inblock_z, 
                                    random = ~ 1 | subject, 
                                    data=dat4Analyses_reducedby_8, 
                                    method="ML",
                                    na.action = na.omit)
# random = ~ 1 | subject/num_block_z, 
summary(fit_dat4Analyses_reducedby_8)
intervals(fit_dat4Analyses_reducedby_8)
reghelper::simple_slopes(fit_dat4Analyses_reducedby_8, 
                         levels = list(num_inblock = 'sstest',
                                       accuracy = c(0, 1)),
                         confint = T)

#####################
# 2) Simple Mean Analysis -- nature of the accuracy main effect
library(emmeans)
fit_dat4Analyses_reducedby_8 <- lme(erp_ern ~ accuracy*observation*num_block*num_inblock, 
                                    random = ~ 1 | subject, 
                                    data=dat4Analyses_reducedby_8, 
                                    method="ML",
                                    na.action = na.omit)
# Obtain estimated marginal means for the main effect of accuracy
accuracy_emm <- emmeans(fit_dat4Analyses_reducedby_8, ~ accuracy)

# Display the estimated marginal means
summary(accuracy_emm)
# accuracy emmean    SE df lower.CL upper.CL
# 0        -2.760 0.352 43   -3.470    -2.05
# 1         0.559 0.331 43   -0.108     1.23

# Perform pairwise comparisons (if there are more than two levels)
pairwise_comparisons <- pairs(accuracy_emm)
summary(pairwise_comparisons)

#####################
# 3) Simple Slopes Analysis -- nature of the three-way interaction 
# 1 -- slopes
# Estimate slopes for num_inblock at the mean of num_block
# setting 'at' for fixed value of num_block and calculating slopes for num_inblock
fit_dat4Analyses_reducedby_8 <- lme(erp_ern_z ~ accuracy*observation*num_block_z*num_inblock_z, 
                                    random = ~ 1 | subject, 
                                    data=dat4Analyses_reducedby_8, 
                                    method="ML",
                                    na.action = na.omit)
trends <- emtrends(fit_dat4Analyses_reducedby_8, ~ accuracy * observation, 
                   "num_inblock_z", 
                   at = list(num_block = mean(dat4Analyses_reducedby_8$num_block_z)))
summary(trends)
tests <- test(trends, null = 0)
summary(tests)

# 2 -- contrasts
# Set up the emmeans for the interaction of interest, specifically looking at slopes
emm_slopes <- emtrends(fit_dat4Analyses_reducedby_8, 
                       specs = pairwise ~ accuracy * observation, 
                       var = "num_inblock_z", 
                       at = list(num_block = mean(dat4Analyses_reducedby_8$num_block_z)))
summary(emm_slopes)

# Pairwise comparisons of slopes
slopes_comparisons <- pairs(emm_slopes)
summary(slopes_comparisons)


########################################## exploratory follow-ups ##########################################
# 1) model for feedback
fit_dat4Analyses_reducedby_8 <- lme(erp_ern_z ~ accuracy*observation*num_inblock_z*go_faster_contrast, 
                                    random = ~ 1 | subject, 
                                    data=dat4Analyses_reducedby_8, 
                                    method="ML",
                                    na.action = na.omit)
summary(fit_dat4Analyses_reducedby_8)


#                                                                Value  Std.Error    DF   t-value p-value
# (Intercept)                                              -0.18636077 0.05310156 11626 -3.509516  0.0005
# accuracy1                                                 0.26217785 0.01252596 11626 20.930764  0.0000
# observation1                                             -0.00739446 0.01246281 11626 -0.593322  0.5530
# num_inblock_z                                            -0.00404637 0.01242494 11626 -0.325666  0.7447
# go_faster_contrast1                                      -0.01127248 0.01381233 11626 -0.816117  0.4144
# accuracy1:observation1                                   -0.00798107 0.01229724 11626 -0.649013  0.5163
# accuracy1:num_inblock_z                                   0.02331275 0.01243533 11626  1.874719  0.0609
# observation1:num_inblock_z                               -0.00161872 0.01242141 11626 -0.130317  0.8963
# accuracy1:go_faster_contrast1                             0.00133833 0.01237201 11626  0.108174  0.9139
# observation1:go_faster_contrast1                          0.01236642 0.01245478 11626  0.992906  0.3208
# num_inblock_z:go_faster_contrast1                         0.00072068 0.01242361 11626  0.058009  0.9537
# accuracy1:observation1:num_inblock_z                      0.02614144 0.01242909 11626  2.103246  0.0355
# accuracy1:observation1:go_faster_contrast1               -0.00421591 0.01229996 11626 -0.342758  0.7318
# accuracy1:num_inblock_z:go_faster_contrast1               0.01046606 0.01243228 11626  0.841846  0.3999
# observation1:num_inblock_z:go_faster_contrast1           -0.00445712 0.01242351 11626 -0.358765  0.7198
# accuracy1:observation1:num_inblock_z:go_faster_contrast1  0.00589841 0.01243227 11626  0.474443  0.6352

# 2) model for rt
dat4Analyses_reducedby_8$log_rt <- log(dat4Analyses_reducedby_8$rt*1000 + 1)
dat4Analyses_reducedby_8$log_rt_z <- scale(dat4Analyses_reducedby_8$log_rt , center = TRUE, scale = TRUE)
fit_dat4Analyses_reducedby_8 <- lme(erp_ern_z ~ accuracy*observation*num_inblock_z*log_rt_z, 
                                    random = ~ 1 | subject, 
                                    data=dat4Analyses_reducedby_8, 
                                    method="ML",
                                    na.action = na.omit)
summary(fit_dat4Analyses_reducedby_8)

# Fixed effects:  erp_ern_z ~ accuracy * observation * num_inblock_z * log_rt_z 
#                                                     Value  Std.Error    DF   t-value p-value
# (Intercept)                                   -0.13765120 0.05122767 14168 -2.687048  0.0072
# accuracy1                                      0.22197192 0.01304775 14168 17.012281  0.0000
# observation1                                  -0.00527196 0.01267664 14168 -0.415880  0.6775
# num_inblock_z                                 -0.00217639 0.01217709 14168 -0.178728  0.8582
# log_rt_z                                       0.01278302 0.01122897 14168  1.138396  0.2550
# accuracy1:observation1                        -0.00799167 0.01257614 14168 -0.635463  0.5251
# accuracy1:num_inblock_z                        0.01378753 0.01218520 14168  1.131498  0.2579
# observation1:num_inblock_z                    -0.01334164 0.01217004 14168 -1.096269  0.2730
# accuracy1:log_rt_z                            -0.08051258 0.00976883 14168 -8.241783  0.0000
# observation1:log_rt_z                          0.00555092 0.00974069 14168  0.569869  0.5688
# num_inblock_z:log_rt_z                         0.01568093 0.00941030 14168  1.666357  0.0957
# accuracy1:observation1:num_inblock_z           0.03297833 0.01217506 14168  2.708678  0.0068
# accuracy1:observation1:log_rt_z               -0.00815124 0.00960489 14168 -0.848655  0.3961
# accuracy1:num_inblock_z:log_rt_z               0.01492230 0.00940321 14168  1.586937  0.1125
# observation1:num_inblock_z:log_rt_z            0.00467014 0.00939858 14168  0.496898  0.6193
# accuracy1:observation1:num_inblock_z:log_rt_z  0.00537482 0.00938916 14168  0.572450  0.5670
