# Title: Functional brain connectivity predictors of prospective substance use initiation and their environmental correlates
# Contact: Omid Kardan omidk@med.umich.edu
# This script uses the rsFC values generated from SUI_abcd_PLS_main_and_supp.m and the SUI_in_ABCD_compile.R files 

# Requires downloaded ABCD data tables: abcd_p_demo.csv, abcd_y_lt.csv, ce_p_nsc.csv, ce_y_nsc.csv, led_l_crime.csv, ce_p_fes.csv,
# ce_y_fes.csv, led_l_no2.csv, led_l_pm25.csv, led_l_particulat.csv, su_y_sui.csv, su_y_mypi.csv, ledl_l_adi.csv, mh_p_fhx.csv

# Also requires abcd_sub_event_list.csv:
# create a file named "abcd_sub_event_list.csv" with two columns. one column containing all ABCD subids repeated in 5 rows per sub (one row for each year) and 2nd column with years labeled
# as 'eventname' with levels: 'baseline_year_1_arm_1', '1_year_follow_up_y_arm_1', '2_year_follow_up_y_arm_1', '3_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1'


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


setwd('working_directory_address/')

df_filled_demog_bcog_Y0 <- read.csv('df_filled_demog_bcog_Y0.csv')
SU_us <- read.csv('abcd_pls_rsFC_SU_multrun_scores_forstudy2.csv')

df_total0 <- merge(SU_us,df_filled_demog_bcog_Y0,by.x = 'subid', by.y = 'src_subject_id')
df_total0z <- data.frame(scale(df_total0[,c('U_y2','U_y0','pFD_y0', 'pFD_y2', 'any_year_SU',
                                            'Male_bin','interview_age', 'Income','HighestEd',
                                            'White' , 'Black' , 'Hispanic'  , 'FamilyConflict',
                                            'famhx_ss_momdad_dg_p' , 'famhx_ss_momdad_alc_p','ADI','NeighCrime','NeighSafety',
                                            'air_pollution_no2','air_pollution_pm25', 'reshist_addr1_pb')]))


df_total0z$site_id_l <- df_total0$site_id_l
df_total0z$used_ids <- df_total0$used_ids
df_total0z$any_year_SU <- df_total0$any_year_SU
df_total0z$V_diff_y0 <- df_total0$V_diff_y0
df_total0z$V_diff_y2 <- df_total0$V_diff_y2
df_total0z$sub_id <- df_total0$subid

df_total2d <- df_total0z[df_total0z$used_ids > 0 & !is.na(df_total0z$used_ids),]
df_total2d$expscore <- df_total2d$U_y2*df_total2d$V_diff_y2 + df_total2d$U_y0*df_total2d$V_diff_y0
df_total2d$meanFD <- df_total2d$pFD_y2 + df_total2d$pFD_y0



ll2d <- lmer(used_ids  ~ 1  + meanFD + Male_bin + interview_age + Income + HighestEd +White + Black + Hispanic+ 
               (1|site_id_l) , data=df_total2d)  
summary(ll2d)  # no difference between the groups in terms of the matched categories including sites
ll2d <- lm(used_ids  ~ 1  + meanFD + Male_bin + interview_age + Income + HighestEd +White + Black + Hispanic 
               , data=df_total2d) 
summary(ll2d)  # no difference between the groups in terms of the matched categories ignoring sites

########################## proj score by types of substances
dat_list <- read.csv('working_directory_address/abcd_sub_event_list.csv') # csv including list of subs repeated 5 rows for each year
su_mypi0 <- read.csv('working_directory_address/su_y_mypi.csv') # load ABCD data table su_y_mypi.csv  
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

su_sui0 <- read.csv('working_directory_address/su_y_sui.csv') # load ABCD data table su_y_sui.csv 

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

SUI_types_df <- sub_use_bySub %>% dplyr::select(src_subject_id, alc, nic, mj, other) %>% right_join(df_total2d[,c("sub_id", "expscore", "used_ids")], by = c('src_subject_id' = "sub_id")) 
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

# ##########  Figure 2 in Study 1 results:

counts <- SUI_types_df %>% ungroup() %>% group_by(Category) %>% dplyr::count()
SUI_types_df$Category <- factor(SUI_types_df$Category, levels = c("Control","Alcohol","Cannabis","Nicotine","Other","Poly"))

p1 <- ggplot(SUI_types_df) + aes(x = Category, y = expscore, color = Category, fill = Category) +
  stat_summary(fun = mean, geom = "col", alpha = .3) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", color = "black", width = .2) + 
  geom_text(aes(label = n, x = Category, y = -0.2), data = counts, color = "black") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1), text=element_text(size=16)) + 
  labs(y = "rsFC Pattern Expression Score")
p1

ggsave(filename = "working_directory_address/RplotFigure2.png", plot = p1, height = 5, width = 7, dpi = 500 )


t.test(SUI_types_df[SUI_types_df$Category == "Control",'expscore'], SUI_types_df[SUI_types_df$Category == "Cannabis",'expscore']) # adj p = 5*0.002
t.test(SUI_types_df[SUI_types_df$Category == "Control",'expscore'], SUI_types_df[SUI_types_df$Category == "Other",'expscore']) # adj p = 4*0.003
t.test(SUI_types_df[SUI_types_df$Category == "Control",'expscore'], SUI_types_df[SUI_types_df$Category == "Poly",'expscore']) # adj p = 3*0.0005
t.test(SUI_types_df[SUI_types_df$Category == "Control",'expscore'], SUI_types_df[SUI_types_df$Category == "Nicotine",'expscore']) # adj p = 2*0.016
t.test(SUI_types_df[SUI_types_df$Category == "Control",'expscore'], SUI_types_df[SUI_types_df$Category == "Alcohol",'expscore']) # adj p = 1*0.256  N.S.

t.test(SUI_types_df[SUI_types_df$Category == "Cannabis",'expscore'], SUI_types_df[SUI_types_df$Category == "Alcohol",'expscore']) # p = 1*0.067  N.S.


################################################## Study 2 regressions

df_total00 <- df_total0z[df_total0z$used_ids <1 & !is.na(df_total0z$used_ids)  &  df_total0z$any_year_SU ==0,] # remainder of participants (not included in study 1)
df_total00$projscore <-  df_total00$U_y2*df_total00$V_diff_y2 + df_total00$U_y0*df_total00$V_diff_y0
df_total00$meanFD <- df_total00$pFD_y2 + df_total00$pFD_y0

# pick the variable of interest (VOI) for the regressions: (comment out the rest)
df_total00$VOI <- df_total00$famhx_ss_momdad_dg_p
df_total00$VOI <- df_total00$famhx_ss_momdad_alc_p
df_total00$VOI <- df_total00$FamilyConflict
df_total00$VOI <- df_total00$air_pollution_pm25  
df_total00$VOI <- df_total00$air_pollution_no2  
df_total00$VOI <- df_total00$reshist_addr1_pb   
df_total00$VOI <- df_total00$ADI
df_total00$VOI <- df_total00$NeighSafety
df_total00$VOI <- df_total00$NeighCrime


# run the regressions
df_total_c <- df_total00[(!is.na(df_total00$projscore) &  !is.na(df_total00$meanFD) 
                          &  !is.na(df_total00$famhx_ss_momdad_dg_p) & !is.na(df_total00$famhx_ss_momdad_alc_p)
                          &  !is.na(df_total00$FamilyConflict) & !is.na(df_total00$air_pollution_pm25)
                          &  !is.na(df_total00$air_pollution_no2) & !is.na(df_total00$reshist_addr1_pb)
                          &  !is.na(df_total00$ADI) & !is.na(df_total00$NeighSafety)
                          &  !is.na(df_total00$NeighCrime) 
                          & !is.na(df_total00$Male_bin) & !is.na(df_total00$interview_age) & !is.na(df_total00$Income)
                          & !is.na(df_total00$HighestEd)),]

# version without race/ethnicity as covariates
ll <- lmer(projscore  ~ 1  + meanFD + Male_bin + interview_age + Income + HighestEd 
           + VOI
           + (1|site_id_l) , data=df_total_c)  
lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}
lm.beta.lmer(ll)
summary(ll)

# version with race/ethnicity as covariates
ll <- lmer(projscore  ~ 1  + meanFD + Male_bin + interview_age + Income + HighestEd + White + Black + Hispanic
           + VOI
           + (1|site_id_l) , data=df_total_c)  

lm.beta.lmer(ll)
summary(ll)

