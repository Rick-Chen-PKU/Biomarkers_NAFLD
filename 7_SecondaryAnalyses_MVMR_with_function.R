#######################################################################
#             Biomarker and NAFLD and its clinical spectrum           #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS

# This script calculate the MVMR results


#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory
setwd("/home/dongzechen/GWAS/04-twosampleMR")   # your_working_directory

library(TwoSampleMR)
library(MendelianRandomization) 
library(tidyverse)
library(MVMR)
library(data.table)
library(writexl)
library(readxl)
library(stringr)
library(LDlinkR)

# Preparation
pre_path =  "/mnt/data1/user/chendongze/Project/Biomarkers_NAFLD_Project/disease_"
# exposure_trait_1 = '2_Ala'    # "2_Ala" "4_HDL" "4_TG2021" "6_Ferr" "6_Iron" "6_TSAT" "7_VD" "8_C4" "8_IL1RA"
# exposure_trait_2= 'BMI'  # BMI, T2D, Hypertension(ukb-b-12493), Hypothyroidism(ukb-b-4226), ALB, ALT, AST, GGT

outcome_trait = 'NAFLD_FinnGen'  
outcome = fread(paste0(pre_path, outcome_trait, '.txt')) %>% mutate(BETA = log(as.numeric(OR)))
outcome_raw = outcome

my_mvmr <- function(exposure_trait_1_list, exposure_trait_2_list){
  for (exposure_trait_1 in exposure_trait_1_list){
    exp_1 = fread(paste0(pre_path, exposure_trait_1, '.txt'), fill=TRUE) %>% mutate(BETA = log(OR))
    exp_1_raw = exp_1
    exp_1_ivs = fread(paste0('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/NewVerIVs/', exposure_trait_1, '.csv')) %>% tidyr::drop_na(N)
    for (exposure_trait_2 in exposure_trait_2_list){
      message(paste0(" ######## ", exposure_trait_1, " to ", exposure_trait_2, " ########"))
      message(paste0(" ######## ", exposure_trait_1, " to ", exposure_trait_2, " ########"))
      exp_2 = fread(paste0(pre_path, exposure_trait_2, '.txt')) %>% mutate(BETA = log(OR))
      exp_2_raw = exp_2
      exp_2_ivs = fread(paste0('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/NewVerIVs/', exposure_trait_2, '.csv'))
      
      ivs = Reduce(union,  list(exp_1_ivs$SNP, exp_2_ivs$SNP)) 
      ivs = data.frame(SNP=ivs)
      ivs_independent = clump_data(ivs, clump_kb = 10000, clump_r2 = 0.01, pop = "EUR")
      # ivs_independent = data.frame(SNP=ivs)
      
      exp_1 = exp_1_raw[exp_1_raw$SNP %in% as.character(ivs_independent$SNP), ]
      exp_1 = subset(exp_1, select=c('SNP', 'BETA', 'SE', 'P'))
      exp_2 = exp_2_raw[exp_2_raw$SNP %in% as.character(ivs_independent$SNP), ]
      exp_2 = subset(exp_2, select=c('SNP', 'BETA', 'SE', 'P'))
      outcome = outcome_raw[outcome_raw$SNP %in% as.character(ivs_independent$SNP), ]
      outcome = subset(outcome, select=c('SNP', 'BETA', 'SE', 'P'))
      
      dat = merge(exp_1, exp_2, by='SNP')
      dat = merge(dat, outcome, by='SNP')
      dat = dat[dat$P > 5e-8, ]   
      Num_ivs = nrow(dat)  
      
      trait_1_specific_num = length(intersect(dat$SNP, exp_1_ivs$SNP))
      trait_2_specific_num = length(intersect(dat$SNP, exp_2_ivs$SNP))
      
      rawdat_mvmr1 = subset(dat, select=c('BETA.x','BETA.y', 'SE.x', 'SE.y', 'BETA', 'SE', 'SNP'))
      names(rawdat_mvmr1) = c('exp1_beta','exp2_beta', 'exp1_se', 'exp2_se', 'outcome_beta', 'outcome_se', 'SNP')
      rawdat_mvmr1 = na.omit(rawdat_mvmr1)
      writexl::write_xlsx(rawdat_mvmr1, paste0('./Project/Biomarkers_NAFLD_Project/mvmr/',exposure_trait_1,'-',exposure_trait_2, '_to_', outcome_trait,'_IVs.xlsx'), col_names = TRUE)  # 保存工具变量
      
      
      F.data <- format_mvmr(BXGs = rawdat_mvmr1[,c(1,2)],
                            BYG = rawdat_mvmr1[,5],
                            seBXGs = rawdat_mvmr1[,c(3,4)],
                            seBYG = rawdat_mvmr1[,6],
                            RSID = rawdat_mvmr1[,7])
      
      sres <- strength_mvmr(r_input = F.data, gencov = 0)
      
      pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
      
      res <- ivw_mvmr(r_input = F.data)
      
      out_df = cbind(data.frame(confounders_mediators=c(exposure_trait_1,exposure_trait_2), outcome = outcome_trait,  
                                Raw_Trait_NumIVs=c(length(exp_1_ivs$SNP), length(exp_2_ivs$SNP)), 
                                Final_Trait_NumIVs=c(trait_1_specific_num, trait_2_specific_num), NumIVs =Num_ivs), as.data.frame(res)[,c(1,2,4)],
                     conditional_F_statistic=c(sres$exposure1, sres$exposure2), Qstat=c(pres$Qstat,pres$Qstat), Qpval=c(pres$Qpval, pres$Qpval))
      
      
      if (file.exists('./Biomarkers_NAFLD_Project_MR_Result_MVMR.xlsx')){
        old_dat = read_excel('./Biomarkers_NAFLD_Project_MR_Result_MVMR.xlsx')
        new_dat = rbind(old_dat, "", out_df)
        writexl::write_xlsx(new_dat,'./Biomarkers_NAFLD_Project_MR_Result_MVMR.xlsx') 
      } else {
        writexl::write_xlsx(out_df,'./Biomarkers_NAFLD_Project_MR_Result_MVMR.xlsx') 
      } 
    
    }
  
    
  }

}

exposure_trait_1_list = c( "2_Ala" ,"4_HDL" ,"4_TG2021", "6_Ferr", "6_Iron" ,"6_TSAT" ,"7_VD", "8_C4" ,"8_IL1RA")
exposure_trait_2_list = c( 'BMI', 'T2D', 'Hypertension', 'Hypothyroidism','ALB', 'ALT', 'AST', 'GGT')
my_mvmr(exposure_trait_1_list, exposure_trait_2_list)


