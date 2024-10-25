# Title: Functional brain connectivity predictors of prospective substance use initiation and their environmental correlates
# Contact: Omid Kardan omidk@med.umich.edu
# This script is the revised version of SUI_in_ABCD_regressions.R for BP:CNNI 
# Same requirements described in SUI_in_ABCD_regressions.R


library(dplyr)
library(ggplot2)
library(tidyverse)
library(ppcor)
library(fastDummies)
library(sjPlot)
library(corrplot)
library(lme4)
library(lmerTest)
library(plyr)
library(mediation)
library(effectsize)


setwd('wd/')

df_filled_demog_bcog_Y0 <- read.csv('wd/df_filled_demog_bcog_Y0.csv')
df_filled_demog_bcog_Y2 <- read.csv('wd/df_filled_demog_bcog_Y2.csv')
df_filled_demog_bcog_Y2$y2age <- df_filled_demog_bcog_Y2$interview_age
df_filled_demog_bcog_Y2$y2pivocab <- df_filled_demog_bcog_Y2$nihtbx_picvocab_uncorrected
df_filled_demog_bcog_Y2$y2read <- df_filled_demog_bcog_Y2$nihtbx_reading_uncorrected
df_filled_demog_bcog_Y2$y2picture <- df_filled_demog_bcog_Y2$nihtbx_picture_uncorrected
df_filled_demog_bcog_Y2$y2flanker <- df_filled_demog_bcog_Y2$nihtbx_flanker_uncorrected
df_filled_demog_bcog_Y2$y2pattern <- df_filled_demog_bcog_Y2$nihtbx_pattern_uncorrected
df_filled_demog_bcog_Y2$y2nbk <- df_filled_demog_bcog_Y2$NBK_acc
df_filled_demog_bcog_Y2$y2INT <- df_filled_demog_bcog_Y2$PF10_INT_lavaan
df_filled_demog_bcog_Y2$y2EXT <- df_filled_demog_bcog_Y2$PF10_EXT_lavaan
df_filled_demog_bcog_Y2$y2PF <- df_filled_demog_bcog_Y2$PF10_lavaan
SU_us <- read.csv('wd/abcd_pls_rsFC_SUscores_and_subs_2024_233_multrun.csv')

df_totalq <- merge(SU_us,df_filled_demog_bcog_Y0,by.x = 'subid', by.y = 'src_subject_id')
df_total0 <- merge(df_totalq,df_filled_demog_bcog_Y2[,c('y2age','y2pivocab',
                                                        'y2read','y2picture','y2flanker','y2pattern',
                                                        'y2nbk','y2INT','y2EXT','y2PF','src_subject_id')],by.x = 'subid', by.y = 'src_subject_id')

mean(df_total0$pFD_y0[df_total0$used_ids == 1 & !is.na(df_total0$used_ids)])
sd(df_total0$pFD_y0[df_total0$used_ids == 1 & !is.na(df_total0$used_ids)])

mean(df_total0$pFD_y2[df_total0$used_ids == 1 & !is.na(df_total0$used_ids)])
sd(df_total0$pFD_y2[df_total0$used_ids == 1& !is.na(df_total0$used_ids)])

mean(df_total0$pFD_y0[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)])
sd(df_total0$pFD_y0[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)])

mean(df_total0$pFD_y2[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)])
sd(df_total0$pFD_y2[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)])



t.test(df_total0$pFD_y0[df_total0$used_ids == 1],df_total0$pFD_y2[df_total0$used_ids == 1],paired = TRUE) # t = 0.76389, df = 227, p-value = 0.4457
t.test(df_total0$pFD_y0[df_total0$used_ids == 2],df_total0$pFD_y2[df_total0$used_ids == 2],paired = TRUE) # t = 1.2262, df = 232, p-value = 0.2214
t.test(df_total0$pFD_y0[df_total0$used_ids == 1],df_total0$pFD_y0[df_total0$used_ids == 2])  # t = -0.8724, df = 458.99, p-value = 0.3834
t.test(df_total0$pFD_y2[df_total0$used_ids == 1],df_total0$pFD_y2[df_total0$used_ids == 2]) # t = -0.60668, df = 458.99, p-value = 0.5444

t.test(df_total0$interview_age[df_total0$used_ids == 1],df_total0$interview_age[df_total0$used_ids == 2])  # t = -0.8724, df = 458.99, p-value = 0.3834
t.test(df_total0$y2age[df_total0$used_ids == 1],df_total0$y2age[df_total0$used_ids == 2]) # t = -0.60668, df = 458.99, p-value = 0.5444


mean(df_total0$interview_age[df_total0$used_ids == 1 & !is.na(df_total0$used_ids)]) 
sd(df_total0$interview_age[df_total0$used_ids == 1 & !is.na(df_total0$used_ids)])

mean(df_total0$y2age[df_total0$used_ids == 1 & !is.na(df_total0$used_ids)])
sd(df_total0$y2age[df_total0$used_ids == 1 & !is.na(df_total0$used_ids)])

mean(df_total0$interview_age[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)])
sd(df_total0$interview_age[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)])

mean(df_total0$y2age[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)])
sd(df_total0$y2age[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)])

mean(df_total0$y2age[df_total0$used_ids == 0 & !is.na(df_total0$used_ids)])
sd(df_total0$y2age[df_total0$used_ids == 0 & !is.na(df_total0$used_ids)])

mean(df_total0$interview_age[df_total0$used_ids == 0 & !is.na(df_total0$used_ids)])
sd(df_total0$interview_age[df_total0$used_ids == 0 & !is.na(df_total0$used_ids)])

income_midpoints <- c(5000, 8500, 14000, 20500, 30000, 42500, 62500, 87500, 150000, 200000) 
mean(df_total0$Income[df_total0$used_ids == 1 & !is.na(df_total0$used_ids)],na.rm = TRUE) # 7.3185
income_midpoints[7]+ .3185*(income_midpoints[8]-income_midpoints[7])
min(df_total0$Income[df_total0$used_ids == 1 & !is.na(df_total0$used_ids)],na.rm=TRUE)
max(df_total0$Income[df_total0$used_ids == 1 & !is.na(df_total0$used_ids)],na.rm=TRUE)

mean(df_total0$Income[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)],na.rm = TRUE) # 7.180
income_midpoints[7]+ .1804*(income_midpoints[8]-income_midpoints[7])
min(df_total0$Income[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)],na.rm=TRUE)
max(df_total0$Income[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)],na.rm=TRUE)

mean(df_total0$Income[df_total0$used_ids == 0 & !is.na(df_total0$used_ids)],na.rm = TRUE) # 7.7153
income_midpoints[7]+ .7153*(income_midpoints[8]-income_midpoints[7])
min(df_total0$Income[df_total0$used_ids == 0 & !is.na(df_total0$used_ids)],na.rm=TRUE)
max(df_total0$Income[df_total0$used_ids == 0 & !is.na(df_total0$used_ids)],na.rm=TRUE)


mean(df_total0$HighestEd[df_total0$used_ids == 1 & !is.na(df_total0$used_ids)],na.rm = TRUE)
sd(df_total0$HighestEd[df_total0$used_ids == 1 & !is.na(df_total0$used_ids)],na.rm = TRUE)

mean(df_total0$HighestEd[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)],na.rm = TRUE)
sd(df_total0$HighestEd[df_total0$used_ids == 2 & !is.na(df_total0$used_ids)],na.rm = TRUE)

mean(df_total0$HighestEd[df_total0$used_ids == 0 & !is.na(df_total0$used_ids)],na.rm = TRUE)
sd(df_total0$HighestEd[df_total0$used_ids == 0 & !is.na(df_total0$used_ids)],na.rm = TRUE)


t.test(df_total0$Male_bin[df_total0$used_ids == 2],df_total0$Male_bin[df_total0$used_ids == 0 ]) # t = -1.3602, df = 265.44, p-value = 0.1749
t.test(df_total0$White[df_total0$used_ids == 2],df_total0$White[df_total0$used_ids == 0 ]) # t = -3.4824, df = 263.43, p-value = 0.0005815
t.test(df_total0$Black[df_total0$used_ids == 2],df_total0$Black[df_total0$used_ids == 0 ]) # t = 2.0481, df = 255.24, p-value = 0.04158
t.test(df_total0$Hispanic[df_total0$used_ids == 2],df_total0$Hispanic[df_total0$used_ids == 0 ]) # t = 1.1776, df = 261.18, p-value = 0.24
t.test(df_total0$Other[df_total0$used_ids == 2],df_total0$Other[df_total0$used_ids == 0 ]) # t = 2.265, df = 254.12, p-value = 0.02436
t.test(df_total0$interview_age[df_total0$used_ids == 2],df_total0$interview_age[df_total0$used_ids == 0 ]) # t = 5.0742, df = 265.9, p-value = 7.315e-07
t.test(df_total0$y2age[df_total0$used_ids == 2],df_total0$y2age[df_total0$used_ids == 0 ]) # t = 3.8811, df = 267.17, p-value = 0.000131
t.test(df_total0$Income[df_total0$used_ids == 2],df_total0$Income[df_total0$used_ids == 0 ]) # t = -3.5357, df = 250.2, p-value = 0.0004845
t.test(df_total0$HighestEd[df_total0$used_ids == 2],df_total0$HighestEd[df_total0$used_ids == 0 ]) # t = -4.5689, df = 251.02, p-value = 7.695e-06

c_cog <- c('nihtbx_picvocab_uncorrected','y2pivocab',
           'nihtbx_reading_uncorrected','y2read','nihtbx_picture_uncorrected','y2picture',
           'nihtbx_flanker_uncorrected','y2flanker','nihtbx_pattern_uncorrected','y2pattern',
           'NBK_acc','y2nbk')
c_pfac <- c('PF10_lavaan','y2PF')
df_total0$meanCog <- rowMeans(df_total0[,c_cog], na.rm = TRUE)
df_total0$meanPfac <- rowMeans(df_total0[,c_pfac], na.rm = TRUE)

t.test(df_total0$meanCog[df_total0$used_ids == 1],df_total0$meanCog[df_total0$used_ids == 2])  # t = 2.1426, df = 438.83, p-value = 0.03269
t.test(df_total0$meanPfac[df_total0$used_ids == 1],df_total0$meanPfac[df_total0$used_ids == 2]) # t = -1.1707, df = 452.97, p-value = 0.2423
cohens_d(df_total0$meanCog[df_total0$used_ids == 1],df_total0$meanCog[df_total0$used_ids == 2]) # cd = 0.20, [0.02, 0.39]
cohens_d(df_total0$meanPfac[df_total0$used_ids == 1],df_total0$meanPfac[df_total0$used_ids == 2]) # cd = -0.11, [-0.29, 0.07]

t.test(df_total0$PF10_lavaan[df_total0$used_ids == 1],df_total0$PF10_lavaan[df_total0$used_ids == 2])  # t = -1.0184, df = 453.82, p-value = 0.309
t.test(df_total0$y2PF[df_total0$used_ids == 1],df_total0$y2PF[df_total0$used_ids == 2]) # t = -1.2124, df = 427.72, p-value = 0.226

t.test(df_total0$PF10_EXT_lavaan[df_total0$used_ids == 1],df_total0$PF10_EXT_lavaan[df_total0$used_ids == 2])  # t = -2.2901, df = 457.82, p-value = 0.02247
t.test(df_total0$y2EXT[df_total0$used_ids == 1],df_total0$y2EXT[df_total0$used_ids == 2]) # t = -3.5516, df = 403.87, p-value = 0.000428

t.test(df_total0$PF10_INT_lavaan[df_total0$used_ids == 1],df_total0$PF10_INT_lavaan[df_total0$used_ids == 2])  # t = 1.3555, df = 448.19, p-value = 0.1759
t.test(df_total0$y2INT[df_total0$used_ids == 1],df_total0$y2INT[df_total0$used_ids == 2]) # t = 1.9502, df = 431.39, p-value = 0.0518

t.test(df_total0$NBK_acc[df_total0$used_ids == 1],df_total0$NBK_acc[df_total0$used_ids == 2])  # t = 0.50356, df = 436.36, p-value = 0.6148
t.test(df_total0$y2nbk[df_total0$used_ids == 1],df_total0$y2nbk[df_total0$used_ids == 2]) # t = 0.076076, df = 432.95, p-value = 0.9394

t.test(df_total0$nihtbx_picvocab_uncorrected[df_total0$used_ids == 1],df_total0$nihtbx_picvocab_uncorrected[df_total0$used_ids == 2])  #t = 0.87193, df = 438.35, p-value = 0.3837
t.test(df_total0$y2pivocab[df_total0$used_ids == 1],df_total0$y2pivocab[df_total0$used_ids == 2]) # t = 3.1484, df = 411.63, p-value = 0.001761

t.test(df_total0$nihtbx_reading_uncorrected[df_total0$used_ids == 1],df_total0$nihtbx_reading_uncorrected[df_total0$used_ids == 2])  #t = 0.32217, df = 438.5, p-value = 0.7475
t.test(df_total0$y2read[df_total0$used_ids == 1],df_total0$y2read[df_total0$used_ids == 2]) # t = 1.2806, df = 419.95, p-value = 0.201

t.test(df_total0$nihtbx_picture_uncorrected[df_total0$used_ids == 1],df_total0$nihtbx_picture_uncorrected[df_total0$used_ids == 2])  #t = 1.2508, df = 437.5, p-value = 0.2117
t.test(df_total0$y2picture[df_total0$used_ids == 1],df_total0$y2picture[df_total0$used_ids == 2]) # t = 0.84952, df = 423.02, p-value = 0.3961

t.test(df_total0$nihtbx_flanker_uncorrected[df_total0$used_ids == 1],df_total0$nihtbx_flanker_uncorrected[df_total0$used_ids == 2])  #t = -0.2168, df = 431.76, p-value = 0.8285
t.test(df_total0$y2flanker[df_total0$used_ids == 1],df_total0$y2flanker[df_total0$used_ids == 2]) # t = 1.1095, df = 425.47, p-value = 0.2678

t.test(df_total0$nihtbx_pattern_uncorrected[df_total0$used_ids == 1],df_total0$nihtbx_pattern_uncorrected[df_total0$used_ids == 2])  #t = 1.512, df = 437.35, p-value = 0.1313
t.test(df_total0$y2pattern[df_total0$used_ids == 1],df_total0$y2pattern[df_total0$used_ids == 2]) # t = 0.81297, df = 424.48, p-value = 0.4167


df_total0z <- data.frame(scale(df_total0[,c('U_y2','FD_y0','U_y0','FD_y2', 'pFD_y0', 'pFD_y2', 'any_year_SU',
                                            'Male_bin','interview_age', 'y2age','Income','HighestEd',
                                            'White' , 'Black' , 'Hispanic', 'Asian','Other','nihtbx_picvocab_uncorrected','y2pivocab',
                                            'nihtbx_reading_uncorrected','y2read','nihtbx_picture_uncorrected','y2picture',
                                            'nihtbx_flanker_uncorrected','y2flanker','nihtbx_pattern_uncorrected','y2pattern',
                                            'NBK_acc','y2nbk','PF10_INT_lavaan','y2INT','PF10_EXT_lavaan','y2EXT','PF10_lavaan','y2PF',
                                            'FamilyConflict', 'famhx_ss_momdad_dg_p' , 'famhx_ss_momdad_alc_p','ADI','NeighCrime','NeighSafety',
                                            'air_pollution_no2','air_pollution_pm25','reshist_addr1_pb'
)]))


df_total0z$site_id_l <- df_total0$site_id_l
df_total0z$used_ids <- df_total0$used_ids
df_total0z$any_year_SU <- df_total0$any_year_SU
df_total0z$sub_id <- df_total0$subid
df_total00 <- df_total0z[df_total0z$used_ids <1 & !is.na(df_total0z$used_ids)  &  df_total0z$any_year_SU ==0,]
df_total00$Uscore <-  (df_total00$U_y2*0.6681 +  df_total00$U_y0* (-0.2165)) - (df_total00$U_y2* (0.2167) +  df_total00$U_y0* (-0.6781))
df_total00$meanFD <- df_total00$pFD_y2 + df_total00$pFD_y0
df_total2d <- df_total0z[df_total0z$used_ids > 0 & !is.na(df_total0z$used_ids),]
df_total2d$Uscore <- (df_total2d$U_y2*0.6681 +  df_total2d$U_y0* (-0.2165)) - (df_total2d$U_y2* (0.2167) +  df_total2d$U_y0* (-0.6781))
df_total2d$meanFD <- df_total2d$pFD_y2 + df_total2d$pFD_y0
df_totals2supp <- df_total0z[df_total0z$used_ids >= 0 & !is.na(df_total0z$used_ids),]
df_totals2supp$Uscore <- (df_totals2supp$U_y2*0.6681 +  df_totals2supp$U_y0* (-0.2165)) - (df_totals2supp$U_y2* (0.2167) +  df_totals2supp$U_y0* (-0.6781))
df_totals2supp$meanFD <- df_totals2supp$pFD_y2 + df_totals2supp$pFD_y0

cor(df_total2d$used_ids,df_total2d$Uscore)

numrows <- (df_total2d[df_total2d$Male_bin > 0 & df_total2d$used_ids == 1,]) # 105
numrows <- (df_total2d[df_total2d$White > 0 & df_total2d$used_ids == 1,]) # 113
numrows <- (df_total2d[df_total2d$Black > 0 & df_total2d$used_ids == 1,]) # 33
numrows <- (df_total2d[df_total2d$Hispanic > 0 & df_total2d$used_ids == 1,])# 48
numrows <- (df_total2d[df_total2d$Asian > 0 & df_total2d$used_ids == 1,]) # 0
numrows <- (df_total2d[df_total2d$Other > 0 & df_total2d$used_ids == 1,]) # 34


numrows <- (df_total2d[df_total2d$Male_bin > 0 & df_total2d$used_ids == 2,]) # 108
numrows <- (df_total2d[df_total2d$White > 0 & df_total2d$used_ids == 2,]) # 116
numrows <- (df_total2d[df_total2d$Black > 0 & df_total2d$used_ids == 2,]) # 34
numrows <- (df_total2d[df_total2d$Hispanic > 0 & df_total2d$used_ids == 2,]) # 49
numrows <- (df_total2d[df_total2d$Asian > 0 & df_total2d$used_ids == 2,]) # 0
numrows <- (df_total2d[df_total2d$Other > 0 & df_total2d$used_ids == 2,]) # 34

df_study2<- df_total00[!is.na(df_total00$interview_age) & !is.na(df_total00$Income) & !is.na(df_total00$HighestEd)
                       & !is.na(df_total00$FamilyConflict) & !is.na(df_total00$famhx_ss_momdad_dg_p)
                       & !is.na(df_total00$famhx_ss_momdad_alc_p)  & !is.na(df_total00$ADI)  & !is.na(df_total00$NeighCrime) & !is.na(df_total00$NeighSafety)
                       & !is.na(df_total00$air_pollution_no2)  & !is.na(df_total00$air_pollution_pm25)  & !is.na(df_total00$reshist_addr1_pb) ,]
numrows <- (df_study2[df_study2$Male_bin > 0 ,]) # 1448
numrows <- (df_study2[df_study2$White > 0 ,]) # 1801
numrows <- (df_study2[df_study2$Black > 0 ,]) # 244
numrows <- (df_study2[df_study2$Hispanic > 0 ,])# 500
numrows <- (df_study2[df_study2$Asian > 0 ,]) # 43
numrows <- (df_study2[df_study2$Other > 0 ,]) # 266

df_study2_supp<- df_totals2supp[!is.na(df_totals2supp$interview_age) & !is.na(df_totals2supp$Income) & !is.na(df_totals2supp$HighestEd)
                                & !is.na(df_totals2supp$FamilyConflict) & !is.na(df_totals2supp$famhx_ss_momdad_dg_p)
                                & !is.na(df_totals2supp$famhx_ss_momdad_alc_p)  & !is.na(df_totals2supp$ADI)  & !is.na(df_totals2supp$NeighCrime) & !is.na(df_totals2supp$NeighSafety)
                                & !is.na(df_totals2supp$air_pollution_no2)  & !is.na(df_totals2supp$air_pollution_pm25)  & !is.na(df_totals2supp$reshist_addr1_pb) ,]
numrows <- (df_study2_supp[df_study2_supp$Male_bin > 0 ,]) # 1676
numrows <- (df_study2_supp[df_study2_supp$White > 0 ,]) # 2051
numrows <- (df_study2_supp[df_study2_supp$Black > 0 ,]) # 302
numrows <- (df_study2_supp[df_study2_supp$Hispanic > 0 ,])# 606
numrows <- (df_study2_supp[df_study2_supp$Asian > 0 ,]) # 43
numrows <- (df_study2_supp[df_study2_supp$Other > 0 ,]) # 330

ll2d <- lmer(used_ids  ~ 1  + pFD_y0 + pFD_y2 + Male_bin + interview_age + Income + HighestEd +White + Black + Hispanic+ 
               (1|site_id_l) , data=df_total2d)  
summary(ll2d)  # no difference between the groups in terms of the matched categories



## main results study 2 regressions
df_study2$VOI <- df_study2$famhx_ss_momdad_dg_p
df_study2$VOI <- df_study2$famhx_ss_momdad_alc_p
df_study2$VOI <- df_study2$FamilyConflict
df_study2$VOI <- df_study2$air_pollution_pm25  #
df_study2$VOI <- df_study2$air_pollution_no2  #
df_study2$VOI <- df_study2$reshist_addr1_pb  # 
df_study2$VOI <- df_study2$ADI
df_study2$VOI <- df_study2$NeighSafety
df_study2$VOI <- df_study2$NeighCrime


ll <- lmer(Uscore  ~ 1  + meanFD + Male_bin + interview_age + Income + HighestEd 
           + VOI
           + (1|site_id_l) , data=df_study2)  
lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}
lm.beta.lmer(ll)
summary(ll)

##
ll <- lmer(Uscore  ~ 1  + meanFD + Male_bin + interview_age + Income + HighestEd + White + Black + Hispanic
           + VOI
           + (1|site_id_l) , data=df_study2)  
lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}
lm.beta.lmer(ll)
summary(ll)

## supp results study 2 not excluding study 1 Ps regressions
df_study2_supp$VOI <- df_study2_supp$famhx_ss_momdad_dg_p
df_study2_supp$VOI <- df_study2_supp$famhx_ss_momdad_alc_p
df_study2_supp$VOI <- df_study2_supp$FamilyConflict
df_study2_supp$VOI <- df_study2_supp$air_pollution_pm25  #
df_study2_supp$VOI <- df_study2_supp$air_pollution_no2  #
df_study2_supp$VOI <- df_study2_supp$reshist_addr1_pb  # 
df_study2_supp$VOI <- df_study2_supp$ADI
df_study2_supp$VOI <- df_study2_supp$NeighSafety
df_study2_supp$VOI <- df_study2_supp$NeighCrime


ll <- lmer(Uscore  ~ 1  + meanFD + Male_bin + interview_age + Income + HighestEd 
           + VOI
           + (1|site_id_l) , data=df_study2_supp)  
lm.beta.lmer(ll)
summary(ll)

ll <- lmer(Uscore  ~ 1  + meanFD + Male_bin + interview_age + Income + HighestEd + White + Black + Hispanic
           + VOI
           + (1|site_id_l) , data=df_study2_supp)  
lm.beta.lmer(ll)
summary(ll)



###########
########################## types of substances
dat_list <- read.csv('wd/abcd_sub_event_list.csv')
su_mypi0 <- read.csv('wd/su_y_mypi.csv') 
su_mypi00 <- su_mypi0 %>% mutate(eventname = case_when(
  eventname == "6_month_follow_up_arm_1" ~ "1_year_follow_up_y_arm_1",
  eventname == "18_month_follow_up_arm_1" ~ "2_year_follow_up_y_arm_1",
  eventname == "30_month_follow_up_arm_1" ~ "3_year_follow_up_y_arm_1",
  eventname == "42_month_follow_up_arm_1" ~ "4_year_follow_up_y_arm_1",
))
su_mypi00[su_mypi00 == 777] = NA

c1 <- c("mypi_alc_drink_used", "mypi_tob_used", "mypi_chew_pst_used", "mypi_mj_used", "mypi_mj_tinc_used", "mypi_mj_synth_used", 
        "mypi_sniff_used", "mypi_pills_used", "mypi_pills_dep_used", "mypi_pr_used", "mypi_cold_used", 'mypi_high_other_used', 'mypi_ecstasy_used', 
        "mypi_hallucinogen_used", "mypi_ketamine_used", "mypi_ghb_used", "mypi_cathinones_used", 'mypi_salvia_used', 'mypi_tobo_used', 'mypi_cigar_used', 
        "mypi_coke_used", "mypi_meth_used", "mypi_heroin_used", "mypi_steroids_used", "mypi_mushroom_used")

su_mypi <- merge(dat_list,su_mypi00[,c('src_subject_id','eventname',c1)],
                 by = c('src_subject_id','eventname'), all.x = TRUE)

c_alc <- c("mypi_alc_drink_used")

c_nic <- c("mypi_tob_used", "mypi_chew_pst_used", 'mypi_tobo_used', 'mypi_cigar_used')

c_mj <- c( "mypi_mj_used", "mypi_mj_tinc_used", "mypi_mj_synth_used")

c_other <- c("mypi_mj_synth_used", "mypi_sniff_used", "mypi_pills_used", "mypi_pills_dep_used", "mypi_pr_used", "mypi_cold_used", 'mypi_high_other_used', 'mypi_ecstasy_used', 
             "mypi_hallucinogen_used", "mypi_ketamine_used", "mypi_ghb_used", "mypi_cathinones_used", 'mypi_salvia_used',  
             "mypi_coke_used", "mypi_meth_used", "mypi_heroin_used", "mypi_steroids_used", "mypi_mushroom_used")


su_mypi$mypi_alc <- su_mypi[,c_alc]
su_mypi$mypi_nic <- rowMeans(su_mypi[,c_nic], na.rm = TRUE)
su_mypi$mypi_mj <- rowMeans(su_mypi[,c_mj], na.rm = TRUE)
su_mypi$mypi_other <- rowMeans(su_mypi[,c_other], na.rm = TRUE)

su_mypiFJ <- su_mypi %>% dplyr::select(src_subject_id, eventname, mypi_alc, mypi_nic, mypi_mj, mypi_other)

su_sui0 <- read.csv('wd/su_y_sui.csv') # 
c2<- c(  'tlfb_alc_use_l','tlfb_cig_use_l','tlfb_ecig_use_l','tlfb_chew_use_l','tlfb_cigar_use_l','tlfb_hookah_use_l','tlfb_pipes_use_l',
         'tlfb_nicotine_use_l','tlfb_mj_use_l','tlfb_blunt_use_l','su_tlfb_vape_mj_fl_use','tlfb_edible_use_l','tlfb_mj_conc_use_l',
         'su_tlfb_vape_mj_oils_use','tlfb_mj_drink_use_l','tlfb_tincture_use_l','tlfb_mj_synth_use_l','tlfb_coc_use_l',
         'tlfb_bsalts_use_l','tlfb_meth_use_l','tlfb_mdma_use_l','tlfb_ket_use_l','tlfb_ghb_use_l','tlfb_opi_use_l','tlfb_lsd_use_l','tlfb_lsd_use_type_l___1',
         'tlfb_lsd_use_type_l___2','tlfb_lsd_use_type_l___3','tlfb_lsd_use_type_l___4','tlfb_lsd_use_type_l___6','tlfb_lsd_use_type_l___7','tlfb_lsd_use_type_l___8',
         'tlfb_lsd_use_type_l___10','tlfb_shrooms_use_l','tlfb_salvia_use_l','tlfb_steroids_use_l','tlfb_bitta_use_l','tlfb_inhalant_use_l',
         'tlfb_inhalant_use_type_l___1','tlfb_inhalant_use_type_l___2','tlfb_inhalant_use_type_l___3','tlfb_inhalant_use_type_l___4','tlfb_inhalant_use_type_l___5',
         'tlfb_inhalant_use_type_l___6','tlfb_inhalant_use_type_l___7','tlfb_inhalant_use_type_l___8','tlfb_amp_use_l','tlfb_tranq_use_l','tlfb_vicodin_use_l',
         'tlfb_cough_use_l','tlfb_other_use_l','su_tlfb_vaped_first_use','tlfb_alc_first_use','tlfb_cig_first_use','tlfb_ecig_first_use','tlfb_chew_first_use',
         'tlfb_cigar_first_use','tlfb_hookah_first_use','tlfb_pipes_first_use','tlfb_nicotine_first_use','tlfb_mj_first_use','tlfb_blunt_first_use',
         'su_tlfb_mj_fl_first_use','tlfb_edible_first_use','su_tlfb_mj_oils_first_use','tlfb_mj_con_first_use','tlfb_mj_drink_first_use','tlfb_tincture_first_use',
         'tlfb_mj_synth_first_use','su_tlfb_cbd_first_use','tlfb_coc_first_use','tlfb_bsalts_first_use','tlfb_meth_first_use','tlfb_mdma_first_use','tlfb_ket_first_use',
         'tlfb_ghb_first_use','tlfb_opi_first_use','tlfb_hall_first_use','tlfb_shrooms_first_use','tlfb_salvia_first_use','tlfb_steroids_first_use',
         'tlfb_bitta_first_use','tlfb_inhalant_first_use','tlfb_amp_first_use','tlfb_tranq_first_use','tlfb_vicodin_first_use','tlfb_cough_first_use',
         'tlfb_other_first_use','tlfb_6mo_lsd_use_type_l___1','tlfb_6mo_lsd_use_type_l___2','tlfb_6mo_lsd_use_type_l___3','tlfb_6mo_lsd_use_type_l___4',
         'tlfb_6mo_lsd_use_type_l___6','tlfb_6mo_lsd_use_type_l___7','tlfb_6mo_lsd_use_type_l___8','tlfb_6mo_lsd_use_type_l___9','tlfb_other_use2_l',
         'xskipout_alc_l','xskipout_tob_l','xskipout_mj_l','xskipout_other_l')


su_sui <- merge(dat_list,su_sui0[,c('src_subject_id','eventname',c2)],
                by = c('src_subject_id','eventname'), all.x = TRUE)

su_sui$anytlfbSU <- rowMeans(su_sui[,c2], na.rm = TRUE)

c2_alc <- c(  'tlfb_alc_use_l', 'tlfb_alc_first_use',  'xskipout_alc_l')

c2_nic <- c(  'tlfb_cig_use_l','tlfb_ecig_use_l','tlfb_chew_use_l','tlfb_cigar_use_l','tlfb_hookah_use_l','tlfb_pipes_use_l',
              'tlfb_nicotine_use_l','su_tlfb_vaped_first_use','tlfb_cig_first_use','tlfb_ecig_first_use','tlfb_chew_first_use',
              'tlfb_cigar_first_use','tlfb_hookah_first_use','tlfb_pipes_first_use','tlfb_nicotine_first_use','xskipout_tob_l')


c2_mj <- c(  'tlfb_mj_use_l','tlfb_blunt_use_l','su_tlfb_vape_mj_fl_use','tlfb_edible_use_l','tlfb_mj_conc_use_l',
             'su_tlfb_vape_mj_oils_use','tlfb_mj_drink_use_l','tlfb_tincture_use_l','tlfb_mj_synth_use_l','tlfb_mj_first_use','tlfb_blunt_first_use',
             'su_tlfb_mj_fl_first_use','tlfb_edible_first_use','su_tlfb_mj_oils_first_use','tlfb_mj_con_first_use','tlfb_mj_drink_first_use','tlfb_tincture_first_use',
             'xskipout_mj_l', "su_tlfb_cbd_first_use")


c2_other <- c(  'tlfb_coc_use_l', 'tlfb_bsalts_use_l','tlfb_meth_use_l','tlfb_mdma_use_l','tlfb_ket_use_l','tlfb_ghb_use_l','tlfb_opi_use_l','tlfb_lsd_use_l','tlfb_lsd_use_type_l___1',
                'tlfb_lsd_use_type_l___2','tlfb_lsd_use_type_l___3','tlfb_lsd_use_type_l___4','tlfb_lsd_use_type_l___6','tlfb_lsd_use_type_l___7','tlfb_lsd_use_type_l___8',
                'tlfb_lsd_use_type_l___10','tlfb_shrooms_use_l','tlfb_salvia_use_l','tlfb_steroids_use_l','tlfb_bitta_use_l','tlfb_inhalant_use_l',
                'tlfb_inhalant_use_type_l___1','tlfb_inhalant_use_type_l___2','tlfb_inhalant_use_type_l___3','tlfb_inhalant_use_type_l___4','tlfb_inhalant_use_type_l___5',
                'tlfb_inhalant_use_type_l___6','tlfb_inhalant_use_type_l___7','tlfb_inhalant_use_type_l___8','tlfb_amp_use_l','tlfb_tranq_use_l','tlfb_vicodin_use_l',
                'tlfb_cough_use_l','tlfb_other_use_l','tlfb_coc_first_use','tlfb_bsalts_first_use','tlfb_meth_first_use','tlfb_mdma_first_use','tlfb_ket_first_use',
                'tlfb_ghb_first_use','tlfb_opi_first_use','tlfb_hall_first_use','tlfb_shrooms_first_use','tlfb_salvia_first_use','tlfb_steroids_first_use',
                'tlfb_bitta_first_use','tlfb_inhalant_first_use','tlfb_amp_first_use','tlfb_tranq_first_use','tlfb_vicodin_first_use','tlfb_cough_first_use',
                'tlfb_other_first_use','tlfb_6mo_lsd_use_type_l___1','tlfb_6mo_lsd_use_type_l___2','tlfb_6mo_lsd_use_type_l___3','tlfb_6mo_lsd_use_type_l___4',
                'tlfb_6mo_lsd_use_type_l___6','tlfb_6mo_lsd_use_type_l___7','tlfb_6mo_lsd_use_type_l___8','tlfb_6mo_lsd_use_type_l___9','tlfb_other_use2_l',
                'xskipout_other_l', "tlfb_mj_synth_first_use")

su_sui$anytlfb_alc <- rowMeans(su_sui[,c2_alc], na.rm = TRUE)
su_sui$anytlfb_nic <- rowMeans(su_sui[,c2_nic], na.rm = TRUE)
su_sui$anytlfb_mj <- rowMeans(su_sui[,c2_mj], na.rm = TRUE)
su_sui$anytlfb_other <- rowMeans(su_sui[,c2_other], na.rm = TRUE)


su_suiFJ <- su_sui %>% dplyr::select(src_subject_id, eventname, anytlfb_alc, anytlfb_nic, anytlfb_mj, anytlfb_other, anytlfbSU )

sub_use <- su_suiFJ %>% full_join(su_mypiFJ, by = c('src_subject_id','eventname')) 

sub_use_bySub <- sub_use %>% filter(eventname %in% c('3_year_follow_up_y_arm_1','4_year_follow_up_y_arm_1')) %>% group_by(src_subject_id) %>% 
  dplyr::summarise(across(anytlfb_alc:mypi_other, ~ mean(.x, na.rm = TRUE))) %>% ungroup()
sub_use_bySub$alc <- rowMeans(sub_use_bySub[,c('anytlfb_alc','mypi_alc')], na.rm = TRUE)
sub_use_bySub$nic <- rowMeans(sub_use_bySub[,c('anytlfb_nic','mypi_nic')], na.rm = TRUE)
sub_use_bySub$mj <- rowMeans(sub_use_bySub[,c('anytlfb_mj','mypi_mj')], na.rm = TRUE)
sub_use_bySub$other <- rowMeans(sub_use_bySub[,c('anytlfb_other','mypi_other')], na.rm = TRUE)

SUI_types_df <- sub_use_bySub %>% dplyr::select(src_subject_id, alc, nic, mj, other) %>% 
  right_join(df_total2d[,c("sub_id", "Uscore", "used_ids")], by = c('src_subject_id' = "sub_id")) 
SUI_types_df[is.na(SUI_types_df)] <- 0
SUI_types_df <- SUI_types_df %>% mutate(Poly = (alc > 0) + (nic > 0) + (mj> 0) + (other > 0))
SUI_types_df <- SUI_types_df %>% mutate(Category = case_when(
  used_ids == 1 ~ "Control",
  Poly > 1  ~ "Poly",
  alc > 0 ~ "Alcohol",
  nic > 0 ~ "Nicotine",
  mj > 0 ~ "Cannabis",
  other > 0 ~ "Other"
))

counts <- SUI_types_df %>% ungroup() %>% group_by(Category) %>% dplyr::count()
SUI_types_df$Category <- factor(SUI_types_df$Category, levels = c("Control","Alcohol","Cannabis","Nicotine","Other","Poly"))


p1 <- ggplot(SUI_types_df) + aes(x = Category, y = Uscore, color = Category, fill = Category) +
  stat_summary(fun = mean, geom = "col", alpha = .3) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", color = "black", width = .2) + 
  geom_text(aes(label = n, x = Category, y = -0.1), data = counts, color = "black") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1), text=element_text(size=16)) + 
  labs(y = "rsFC pattern weighted sum score") +coord_cartesian(ylim = c(-.12,.65))
p1
ggsave(filename = "wd/RplotFigure3.png", plot = p1, height = 6, width = 4.5, dpi = 500 )


t.test(SUI_types_df[SUI_types_df$Category == "Control",'Uscore'], SUI_types_df[SUI_types_df$Category == "Other",'Uscore']) # adj p = 5*0.00003
t.test(SUI_types_df[SUI_types_df$Category == "Control",'Uscore'], SUI_types_df[SUI_types_df$Category == "Cannabis",'Uscore']) # adj p = 4*0.00328
t.test(SUI_types_df[SUI_types_df$Category == "Control",'Uscore'], SUI_types_df[SUI_types_df$Category == "Poly",'Uscore']) # adj p = 3*0.0001
t.test(SUI_types_df[SUI_types_df$Category == "Control",'Uscore'], SUI_types_df[SUI_types_df$Category == "Nicotine",'Uscore']) # adj p = 2*0.0005
t.test(SUI_types_df[SUI_types_df$Category == "Control",'Uscore'], SUI_types_df[SUI_types_df$Category == "Alcohol",'Uscore']) # adj p = 1*0.039


SUI_types_df_covs <- SUI_types_df  %>% 
  right_join(df_total2d, by = c('src_subject_id' = "sub_id")) 



c_cog <- c('nihtbx_picvocab_uncorrected','y2pivocab',
           'nihtbx_reading_uncorrected','y2read','nihtbx_picture_uncorrected','y2picture',
           'nihtbx_flanker_uncorrected','y2flanker','nihtbx_pattern_uncorrected','y2pattern',
           'NBK_acc','y2nbk')
c_pfac <- c('PF10_lavaan','y2PF')
SUI_types_df_covs$meanCog <- rowMeans(SUI_types_df_covs[,c_cog], na.rm = TRUE)
SUI_types_df_covs$meanPfac <- rowMeans(SUI_types_df_covs[,c_pfac], na.rm = TRUE)
SUI_types_df_covsc <- SUI_types_df_covs[!is.na(SUI_types_df_covs$interview_age) & !is.na(SUI_types_df_covs$Income) & !is.na(SUI_types_df_covs$HighestEd)
                                        & !is.na(SUI_types_df_covs$meanCog) & !is.na(SUI_types_df_covs$meanPfac),]
countsc <- SUI_types_df_covsc %>% ungroup() %>% group_by(Category) %>% dplyr::count()
ll1 <- lmer(Uscore.x  ~ 1  + meanFD + Male_bin + interview_age + Income + HighestEd + White + Black + Hispanic + meanCog + meanPfac
            + (1|site_id_l) , data=SUI_types_df_covsc)  
lm.beta.lmer(ll1)
summary(ll1)
SUI_types_df_covsc$resid_U <- resid(ll1) + mean(SUI_types_df_covsc$Uscore.x)
p2 <- ggplot(SUI_types_df_covsc) + aes(x = Category, y = resid_U, color = Category, fill = Category) +
  stat_summary(fun = mean, geom = "col", alpha = .3) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", color = "black", width = .2) + 
  geom_text(aes(label = n, x = Category, y = -0.1), data = countsc, color = "black") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1), text=element_text(size=16)) + 
  labs(y = "rsFC pattern weighted sum score (adjusted)") +coord_cartesian(ylim = c(-.12,.65))
p2
ggsave(filename = "wd/RplotFigure3_partial.png", plot = p2, height = 6, width = 4.5, dpi = 500 )


t.test(SUI_types_df_covsc[SUI_types_df_covsc$Category == "Control",'resid_U'], SUI_types_df_covsc[SUI_types_df_covsc$Category == "Other",'resid_U']) # p = 0.040
t.test(SUI_types_df_covsc[SUI_types_df_covsc$Category == "Control",'resid_U'], SUI_types_df_covsc[SUI_types_df_covsc$Category == "Cannabis",'resid_U']) # p = 0.00357
t.test(SUI_types_df_covsc[SUI_types_df_covsc$Category == "Control",'resid_U'], SUI_types_df_covsc[SUI_types_df_covsc$Category == "Poly",'resid_U']) # p = 0.0004
t.test(SUI_types_df_covsc[SUI_types_df_covsc$Category == "Control",'resid_U'], SUI_types_df_covsc[SUI_types_df_covsc$Category == "Nicotine",'resid_U']) # p = 0.0597
