# Title: Functional brain connectivity predictors of prospective substance use initiation and their environmental correlates
# Contact: Omid Kardan omidk@med.umich.edu
# This script compiles all non-brain variables from the ABCD curated files (R5.0)

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
library(plyr)

setwd("working_directory_address/") 
dat_list <- read.csv('working_directory_address/abcd_sub_event_list.csv') # a csv file containing all ABCD subids repeated in 4 rows per sub (one row for each year with years labeled as 'eventname')
##################### demog and race/ethnicity #############################
datafordemo0 <- read.csv('working_directory_address/abcd_p_demo.csv')
datafordemo <- merge(dat_list,datafordemo0[,c('src_subject_id','eventname','demo_comb_income_v2_l',
                                              'demo_prnt_ed_v2_l','demo_prtnr_ed_v2_l',
                                              'race_ethnicity','demo_sex_v2')],
                     by = c('src_subject_id','eventname'), all.x = TRUE)
datafordemo$Income <- datafordemo$demo_comb_income_v2_l
datafordemo$Income[datafordemo$Income==999] <- NA
datafordemo$Income[datafordemo$Income==777] <- NA
datafordemo$Income_cat = factor( datafordemo$Income, levels= 1:10, 
                                 labels = c("5000", "8500", "14000", "20500", "30000",
                                            "42500", "62500", "87500", "150000", "200000") )
datafordemo$HighestEdParent <- datafordemo$demo_prnt_ed_v2_l 
datafordemo$HighestEdParent <- as.numeric(datafordemo$HighestEdParent)
datafordemo$HighestEdParent[datafordemo$HighestEdParent == 999] <- NA
datafordemo$HighestEdParent[datafordemo$HighestEdParent == 777] <- NA
datafordemo$HighestEdPartner <- datafordemo$demo_prtnr_ed_v2_l
datafordemo$HighestEdPartner <- as.numeric(datafordemo$HighestEdPartner)
datafordemo$HighestEdPartner[datafordemo$HighestEdPartner == 999] <- NA
datafordemo$HighestEdPartner[datafordemo$HighestEdPartner == 777] <- NA
#retain the highest education out of the two parents/partners
datafordemo$HighestEd <- pmax(datafordemo$HighestEdParent, datafordemo$HighestEdPartner, na.rm =TRUE)
datafordemo$Male_bin = ifelse(datafordemo$demo_sex_v2 == 1 | datafordemo$demo_sex_v2 == 3,1,0) # 3 intsex_male and no intsex_fem
#dummy code race
datafordemo <- datafordemo %>% 
  dummy_cols(ignore_na = TRUE, select_columns = c("race_ethnicity"))
colnames(datafordemo)[14] = c("White") 
colnames(datafordemo)[15] = c("Black")
colnames(datafordemo)[16] = c("Hispanic")
colnames(datafordemo)[17] = c("Asian")
colnames(datafordemo)[18] = c("Other")

datafordem <- datafordemo[,c('src_subject_id','eventname','Income','HighestEd','Male_bin',
                             'White','Black','Hispanic','Asian','Other')]

##################### age and site and relatives #############################
dataforage0 <- read.csv('working_directory_address/abcd_y_lt.csv')
dataforage <- merge(dat_list,dataforage0[,c('src_subject_id','eventname','site_id_l',
                                            'interview_age','rel_family_id')],
                    by = c('src_subject_id','eventname'), all.x = TRUE)

##################### neigh safety and family conflict #############################
dataforsafe0 <- read.csv('working_directory_address/ce_p_nsc.csv')
dataforsafep <- merge(dat_list,dataforsafe0[,c('src_subject_id','eventname','nsc_p_ss_mean_3_items')],
                      by = c('src_subject_id','eventname'), all.x = TRUE)

dataforsafe0 <- read.csv('working_directory_address/ce_y_nsc.csv')
dataforsafey <- merge(dat_list,dataforsafe0[,c('src_subject_id','eventname','neighborhood_crime_y')],
                      by = c('src_subject_id','eventname'), all.x = TRUE)

dataforcrime0 <- read.csv('working_directory_address/led_l_crime.csv')
dataforcrime <- merge(dat_list,dataforcrime0[,c('src_subject_id','eventname','reshist_addr1_p1vlnt',
                                                'reshist_addr1_drgsale')],
                      by = c('src_subject_id','eventname'), all.x = TRUE)

dataforConf0 <- read.csv('working_directory_address/ce_p_fes.csv')
dataforConfp <- merge(dat_list,dataforConf0[,c('src_subject_id','eventname','fes_p_ss_fc_pr')],
                      by = c('src_subject_id','eventname'), all.x = TRUE)

dataforConf0 <- read.csv('working_directory_address/ce_y_fes.csv')
dataforConfy <- merge(dat_list,dataforConf0[,c('src_subject_id','eventname','fes_y_ss_fc_pr')],
                      by = c('src_subject_id','eventname'), all.x = TRUE)

df_cfc0 <- dataforsafep %>% full_join(dataforcrime, by = c('src_subject_id','eventname')) %>% 
  full_join(dataforConfp, by = c('src_subject_id','eventname')) %>% full_join(dataforsafey, by = c('src_subject_id','eventname'))%>%
  full_join(dataforConfy, by = c('src_subject_id','eventname'))

#youth reports of neigh safety + adult reports of neigh safety
df_cfc0$NeighSafety     <- (df_cfc0$neighborhood_crime_y + df_cfc0$nsc_p_ss_mean_3_items)
#log + 1 of adultviolentcrime + drugsale 
df_cfc0$NeighCrime      <- (log1p(df_cfc0$reshist_addr1_p1vlnt) + log1p(df_cfc0$reshist_addr1_drgsale))
#youth reports of fam conflict + adult reports of fam conflict
df_cfc0$FamilyConflict  <- (df_cfc0$fes_y_ss_fc_pr + df_cfc0$fes_p_ss_fc_pr)
df_cfc <- df_cfc0[,c('src_subject_id','eventname','NeighSafety','NeighCrime','FamilyConflict')]

##################### environment phys #############################

dataforNO20 <- read.csv('working_directory_address/led_l_no2.csv')
dataforNO2 <- merge(dat_list,dataforNO20[,c('src_subject_id','eventname','reshist_addr1_no2_2016_aavg')],
                    by = c('src_subject_id','eventname'), all.x = TRUE)


dataforPM250 <- read.csv('working_directory_address/led_l_pm25.csv')
dataforPM25 <- merge(dat_list,dataforPM250[,c('src_subject_id','eventname','reshist_addr1_pm252016aa')],
                     by = c('src_subject_id','eventname'), all.x = TRUE)

dataforlead01 <- read.csv('working_directory_address/led_l_particulat.csv')
dataforlead1 <- merge(dat_list,dataforlead01[,c('src_subject_id','eventname','reshist_addr1_pb')],
                      by = c('src_subject_id','eventname'), all.x = TRUE)


############SU Use #######################################################
su_mypi0 <- read.csv('working_directory_address/su_y_mypi.csv') 
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

su_mypi$mypiSU <- rowMeans(su_mypi[,c1], na.rm = TRUE)

su_mypiFJ <- su_mypi %>% dplyr::select(src_subject_id, eventname, mypiSU)

summary(su_mypiFJ$mypiSU)


su_sui0 <- read.csv('working_directory_address/su_y_sui.csv') # 
c2<- c(
  'tlfb_alc_use_l','tlfb_cig_use_l','tlfb_ecig_use_l','tlfb_chew_use_l','tlfb_cigar_use_l','tlfb_hookah_use_l','tlfb_pipes_use_l',
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


su_suiFJ <- su_sui %>% dplyr::select(src_subject_id, eventname, anytlfbSU)
summary(su_suiFJ$anytlfbSU)

sub_use <- su_suiFJ %>% full_join(su_mypiFJ, by = c('src_subject_id','eventname')) 

##################### environment ADI #############################
adi0 <- read.csv('working_directory_address/led_l_adi.csv')
library(missMDA)

adi0$ncarlog <- log1p(as.numeric(adi0$reshist_addr1_adi_ncar))
adi0$ntellog <- log1p(as.numeric(adi0$reshist_addr1_adi_ntel))
adi0$nplumblog <- log1p(as.numeric(adi0$reshist_addr1_adi_nplumb))
adi0$crowdlog <- log1p(as.numeric(adi0$reshist_addr1_adi_crowd))

dat <- adi0 %>% dplyr::select(c("reshist_addr1_adi_edu_l" ,"reshist_addr1_adi_edu_h", 
                                "reshist_addr1_adi_work_c","reshist_addr1_adi_income", "reshist_addr1_adi_in_dis",
                                "reshist_addr1_adi_home_v", "reshist_addr1_adi_rent","reshist_addr1_adi_mortg",
                                "reshist_addr1_adi_home_o", "reshist_addr1_adi_unemp", "reshist_addr1_adi_pov",
                                "reshist_addr1_adi_b138", "reshist_addr1_adi_sp" , "ncarlog", 
                                "ntellog", "nplumblog", "crowdlog"))


#nb <- estim_ncpPCA(dat,method.cv = "Kfold", verbose = FALSE)  # nb$ncp = 5
comp <- imputePCA(dat, ncp = 5,method="Regularized")  # nb$ncp = 5
tt<- prcomp(comp$completeObs, center = TRUE, scale = TRUE)
summary(tt)
adi0$ADI_pc1 <- tt$x[,1]

adi <- merge(dat_list,adi0[,c('src_subject_id','eventname','reshist_addr1_adi_wsum','ADI_pc1')],
             by = c('src_subject_id','eventname'), all.x = TRUE)

##################### asq parental use #############################
famhx0 <- read.csv('working_directory_address/mh_p_fhx.csv')
famhx <- merge(dat_list,famhx0[,c('src_subject_id','eventname','famhx_ss_momdad_dg_p','famhx_ss_momdad_alc_p')],
               by = c('src_subject_id','eventname'), all.x = TRUE)


########## combine all data (no cog)
df_total_nocog <- adi %>% 
  full_join(sub_use, by = c('src_subject_id','eventname'))%>%
  full_join(df_cfc, by = c('src_subject_id','eventname'))%>%
  full_join(dataforage, by = c('src_subject_id','eventname'))%>%
  full_join(dataforlead1, by = c('src_subject_id','eventname'))%>%
  full_join(dataforNO2, by = c('src_subject_id','eventname'))%>%
  full_join(datafordem, by = c('src_subject_id','eventname'))%>%
  full_join(dataforPM25, by = c('src_subject_id','eventname'))%>%
  full_join(famhx, by = c('src_subject_id','eventname'))
write.csv(df_total_nocog,'working_directory_address/df_total_nonbrain2024.csv')  


##### generate the non-brain csv files used in SUI_abcd_PLS_main_and_supp.m script
df_total_nocog <- read.csv('working_directory_address/df_total_nonbrain2024.csv') 



dy <- df_total_nocog %>% mutate(Income=as.numeric(Income),HighestEd=as.numeric(HighestEd), Male_bin=as.numeric(Male_bin), White=as.numeric(White),
                                Black=as.numeric(Black), Hispanic=as.numeric(Hispanic), Asian=as.numeric(Asian), Other=as.numeric(Other),
                                lead_1 = as.numeric(reshist_addr1_pb), air_pollution_no2 = as.numeric(reshist_addr1_no2_2016_aavg),air_pollution_pm25 = as.numeric(reshist_addr1_pm252016aa),
                                NeighSafety = as.numeric(NeighSafety),NeighCrime = as.numeric(NeighCrime),FamilyConflict = as.numeric(FamilyConflict),
                                ADI = as.numeric(reshist_addr1_adi_wsum),ADI_pc1 = as.numeric(ADI_pc1))

dy <- dy %>% group_by(src_subject_id) %>% mutate_at(vars(Income, HighestEd, Male_bin, White, Black, Hispanic, Asian, Other,
                                                         lead_1,air_pollution_no2,air_pollution_pm25,NeighSafety,NeighCrime,FamilyConflict,ADI,ADI_pc1),
                                                    ~replace_na(.,mean(.,na.rm=TRUE))) # average the SES measures if they are available for multiple years
write.csv(dy[dy$eventname == 'baseline_year_1_arm_1',],'working_directory_address/df_filled_demog_bcog_Y0.csv') 
write.csv(dy[dy$eventname == '1_year_follow_up_y_arm_1',],'working_directory_address/df_filled_demog_bcog_Y1.csv') 
write.csv(dy[dy$eventname == '2_year_follow_up_y_arm_1',],'working_directory_address/df_filled_demog_bcog_Y2.csv') 
write.csv(dy[dy$eventname == '3_year_follow_up_y_arm_1',],'working_directory_address/df_filled_demog_bcog_Y3.csv') 
write.csv(dy[dy$eventname == '4_year_follow_up_y_arm_1',],'working_directory_address/df_filled_demog_bcog_Y4.csv')
