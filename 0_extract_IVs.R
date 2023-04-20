#######################################################################
#             Biomarker and NAFLD and its clinical spectrum           #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS


#---------------------------------------------------------------------#
# 						Usage Description							  #
#---------------------------------------------------------------------#
# This script is used to extract the IVs.
# Three different approaches were taken to extract IVs, depending on the different data sources and the likelihood of obtaining complete summary statistics. See appendix table for details of data sources, availability of full summary statistic, etc. for each phenotype.
# Method 1 applies to exposed phenotypes for which full summary statistic data is available.
# Method 2 applies to obtaining IVs from the IEU OpenGWAS website.
# Method 3 For manual extraction from the literature when the full summary statistic is neither available nor available from the IEU OpenGWAS website.
# Note: The following code may not work directly, so make simple changes to suit your environment


#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls())   #  Remove any existing objects in R 

# Set working directory 
setwd('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/')   # your_working_directory


#---------------------------------------------------------------------#
#                          Extract  the  IVs                          #----
#---------------------------------------------------------------------#

#------------------------------------------------------------#
#              			Method 1                 			 #----
#------------------------------------------------------------#
library(data.table)
library(stringr)
library(TwoSampleMR)
library(dplyr)

n = 8961   # sample size
# clump parameter 
disease_data <- fread("./disease_Hypothyroidism.txt", na.strings="")  
disease_data = disease_data %>% mutate(SNP=MarkerName, A1=toupper(Allele1), A2=toupper(Allele2), 
									OR=exp(as.numeric(Effect)), SE=as.numeric(StdErr), P=as.numeric(`P-value`), EAF=Freq1, N=n)
disease_data = disease_data %>% select(SNP,A1,A2,OR,SE,P,EAF,N) 
disease_data1 = disease_data
write.table(disease_data, './disease_PA.txt',col.names=TRUE, row.names = FALSE, sep=" ",quote=FALSE)

disease_data = disease_data1
disease_data = disease_data %>% filter(P<5E-6)
disease_data = disease_data %>% mutate(BETA=log(disease_data$OR))
disease_data$Trait = "Hypothyroidism"

exp_dat <- format_data( temp,
						type = "exposure",
						eaf_col = "EAF",
						phenotype_col = "Trait",
						snp_col = "SNP",
						beta_col = "BETA",
						se_col = "SE",
						effect_allele_col = "A1",
						other_allele_col = "A2",
						pval_col = "P",
						samplesize_col = "N"
						)
kb = 10000
r2 = 0.001
p1 = 5E-8
p2 = 5E-3
exp_dat = clump_data(exp_dat, clump_kb = kb, clump_r2 = r2, clump_p1 = p1, clump_p2 = p2, pop = "EUR"); nrow(exp_dat)

if (!is.na(exp_dat$eaf.exposure[1])){
	exp_dat = subset(exp_dat, select=c("SNP", "effect_allele.exposure","other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "exposure", "eaf.exposure"))
	names(exp_dat) = c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'N', 'Trait', 'EAF')
} else {
	exp_dat = subset(exp_dat, select=c("SNP", "effect_allele.exposure","other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "exposure"))
	names(exp_dat) = c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'N', 'Trait')
}
write.csv(exp_dat, paste0('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/TG2021.csv'), row.names = FALSE, quote = FALSE)

#------------------------------------------------------------#
#              			Method 2                 			 #----
#------------------------------------------------------------#
kb = 10000
r2 = 0.001
p1 = 5E-8
p2 = 5E-3
id = 'ebi-a-GCST90012048'  # come from IEU OpenGWAS website
exp_dat = TwoSampleMR::extract_instruments(id, p1 = p1, clump = TRUE, p2 = p2, kb=kb, r2=r2); nrow(exp_dat)   # GWAS ID in IEU OpenGWAS website

if (!is.na(exp_dat$eaf.exposure[1])){
	exp_dat = subset(exp_dat, select=c("SNP", "effect_allele.exposure","other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "exposure", "eaf.exposure"))
	names(exp_dat) = c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'N', 'Trait', 'EAF')
} else {
	exp_dat = subset(exp_dat, select=c("SNP", "effect_allele.exposure","other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "exposure"))
	names(exp_dat) = c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'N', 'Trait')
}
write.csv(exp_dat, paste0('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/TM.csv'), row.names = FALSE, quote = FALSE)


#------------------------------------------------------------#
#              			Method 3                 			 #----
#------------------------------------------------------------#
# Extract directed from the source paper. Ensure P value < 5e-8, kb = 10000 and r2 < 0.001
library(data.table)
library(stringr)
library(TwoSampleMR)
my_trait = "HbA1c"
disease_data <- fread(paste0('./', my_trait, '.csv'), na.strings="")  
disease_data$SNP = gsub(" ", "", disease_data$SNP)   # remove blank in SNP column
exp_dat <- format_data( disease_data,
						type = "exposure",
						eaf_col = "EAF",
						phenotype_col = "Trait",
						snp_col = "SNP",
						beta_col = "BETA",
						se_col = "SE",
						effect_allele_col = "A1",
						other_allele_col = "A2",
						pval_col = "P",
						samplesize_col = "N"
						)

kb = 10000
r2 = 0.001
p1 = 5E-8
P2 = 5E-3
exp_dat = clump_data(exp_dat, clump_kb = kb, clump_r2 = r2, clump_p1 = p1, clump_p2 = p2, pop = "EUR")

if (!is.na(exp_dat$eaf.exposure[1])){
	exp_dat = subset(exp_dat, select=c("SNP", "effect_allele.exposure","other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "exposure", "eaf.exposure"))
	names(exp_dat) = c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'N', 'Trait', 'EAF')
} else {
	exp_dat = subset(exp_dat, select=c("SNP", "effect_allele.exposure","other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "exposure"))
	names(exp_dat) = c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'N', 'Trait')
}
write.csv(exp_dat, paste0('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/', my_trait, '.csv'), row.names = FALSE, quote = FALSE)


