#######################################################################
#             Biomarker and NAFLD and its clinical spectrum           #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS


# This script conduct Harmonisation for main UKB, produce IV Output in two-sample MR format

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory
setwd("/home/dongzechen/GWAS/04-twosampleMR/")  # your_working_directory

if (!require("pacman")){install.packages("pacman")}
pacman::p_load("TwoSampleMR",  "mr.raps", "MendelianRandomization", "MRPRESSO", "data.table", 
				"stringr", "readxl", "writexl", "stringr", "tidyverse")

#---------------------------------------------------------------------#
#                            Harmonisation                            #----
#---------------------------------------------------------------------#
Project = 'Biomarkers_NAFLD_Project'
path = paste0("/mnt/data1/user/chendongze/Project/", Project, "/disease_")


Harmonisation = function(exposure_list, outcome_trait){ 
	# Reading outcome data
	disease_path = paste0(path, outcome_trait,c('.txt'))
	print("Reading outcome data")
	outcome = fread(disease_path) %>% mutate(BETA = log(as.numeric(OR)), Trait = outcome_trait)
	outcome1 = outcome
	
	for (exposure_trait in exposure_list){
		# Reading exposure data
		print(paste0('------Exposure:', exposure_trait, '------'))
		exposure = fread(paste0("./Manual_IVs/", Project, "/NewVerIVs/", exposure_trait, ".csv")) 
		if ('EAF' %in% colnames(exposure)){
			exp_dat = read_exposure_data(
				filename = paste0("./Manual_IVs/", Project, "/NewVerIVs/", exposure_trait, ".csv"),
				sep = ",",
				phenotype_col = "Trait",
				snp_col = "SNP",
				eaf_col = "EAF",
				beta_col = "BETA",
				se_col = "SE",
				effect_allele_col = "A1",
				other_allele_col = "A2",
				pval_col = "P",
				samplesize_col = "N")    
		} else {
			exp_dat = read_exposure_data(
				filename = paste0("./Manual_IVs/", Project, "/NewVerIVs/", exposure_trait, ".csv"),
				sep = ",",
				phenotype_col = "Trait",
				snp_col = "SNP",
				beta_col = "BETA",
				se_col = "SE",
				effect_allele_col = "A1",
				other_allele_col = "A2",
				pval_col = "P",
				samplesize_col = "N")
		}

		
		# Search the LDlink database; Information extract from the former section "2_find_LD_proxies_for_outcome"
		if (outcome_trait=='NAFLD_Anstee'){
			proxy_snp = proxy_snp_anstee
		} else if (outcome_trait=='NAFLD_FinnGen'){
			proxy_snp = proxy_snp_finn
		} else if (outcome_trait=='NAFLD_UKB'){
			proxy_snp = proxy_snp_ukb
		}
		
		intersect_snp = intersect(exp_dat$SNP, unname(proxy_snp))
		if (length(intersect_snp) != 0){
			idx = which(proxy_snp %in% intersect_snp)
			replace_snp = names(proxy_snp)[idx]; names(replace_snp) = intersect_snp
			exp_dat$SNP = str_replace_all(exp_dat$SNP, replace_snp)
			names(intersect_snp) = unname(replace_snp)
		}
		
		outcome = outcome1    # Important
		if ('EAF' %in% colnames(outcome)){
			outcome = subset(outcome, select=c('SNP','A1', 'A2', 'BETA', 'SE', 'P', 'N', 'EAF', 'Trait'))
			print("Writing the outcome to csv")
			outcome = outcome[outcome$SNP %in% exp_dat$SNP, ]
			write.csv(outcome,"./outcome.csv",row.names = FALSE,quote = FALSE)   # Takes a long
			print("Clumping the outcome, take some time")
			outcome_dat <- read_outcome_data(
				filename = "./outcome.csv",
				sep = ",",
				phenotype_col= "Trait",
				snp_col = "SNP",
				eaf_col = "EAF",
				beta_col = "BETA",
				se_col = "SE",
				effect_allele_col = "A1",
				other_allele_col = "A2",
				pval_col = "P",
				samplesize_col = "N")
		} else {
			outcome = subset(outcome, select=c('SNP','A1', 'A2', 'BETA', 'SE', 'P', 'N', 'Trait'))
			print("Writing the outcome to csv")
			outcome = outcome[outcome$SNP %in% exp_dat$SNP, ]
			write.csv(outcome,"./outcome.csv",row.names = FALSE, quote = FALSE)   
			print("Clumping the outcome, take some time")
			outcome_dat <- read_outcome_data(
				filename = "./outcome.csv",
				sep = ",",
				phenotype_col= "Trait",
				snp_col = "SNP",
				beta_col = "BETA",
				se_col = "SE",
				effect_allele_col = "A1",
				other_allele_col = "A2",
				pval_col = "P",
				samplesize_col = "N")
		}
		
		dat = TwoSampleMR::harmonise_data(exp_dat, outcome_dat, action=2); nrow(dat)  # defalut action=2 would remove more SNPs, Remove SNPs for being palindromic with intermediate allele frequencies
		if (length(intersect_snp) != 0){
			dat$SNP = str_replace_all(dat$SNP, intersect_snp)
		}
		# print(paste0('---###palindromic remove SNPs: ', paste(dat[dat$ambiguous == 'TRUE', ]$SNP, collapse=", "), '; total of', length(dat[dat$ambiguous == 'TRUE', ]$SNP)))
		dat = dat[dat$ambiguous == 'FALSE', ]; nrow(dat)  # remove ambiguous SNPs
		# dat$samplesize.outcome = 213431   # sometimes for IEU OPEN GWAS dat
		dat = TwoSampleMR::steiger_filtering(dat) # steiger_filtering, must contain sample_size column
		# print(paste0('---###steiger_filtering remove SNPs: ', paste(dat[dat$steiger_dir == 'FALSE', ]$SNP, collapse=", "), '; total of', length(dat[dat$steiger_dir == 'FALSE', ]$SNP)))
		dat = dat[dat$steiger_dir == 'TRUE', ]; nrow(dat) # remove SNPs didn't pass steiger filtering test
		dat = dat[!is.na(dat$SNP), ]  # remove NA snp
		write.table(dat, paste0('./Project/',Project, '/', exposure_trait, '_to_', outcome_trait, '_instruments.csv'),
			col.names=TRUE, row.names = FALSE, sep=",", quote=FALSE)
	}		
}

exposure_list = list.files('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/NewVerIVs')
exposure_list = gsub(".csv", "", exposure_list)

proxy_snp_1 = read_excel('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/proxy_snp_anstee.xlsx')
proxy_snp_anstee = proxy_snp_1$need_proxy; names(proxy_snp_anstee) = proxy_snp_1$proxy_snp
proxy_snp_2 = read_excel('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/proxy_snp_ukb.xlsx')
proxy_snp_ukb = proxy_snp_2$need_proxy; names(proxy_snp_ukb) = proxy_snp_2$proxy_snp
proxy_snp_3 = read_excel('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/proxy_snp_finn.xlsx')
proxy_snp_finn = proxy_snp_3$need_proxy; names(proxy_snp_finn) = proxy_snp_3$proxy_snp

outcome_list = c('NAFLD_Anstee', 'NAFLD_FinnGen', 'NAFLD_UKB')

for (outcome_trait in outcome_list){
	print(paste0('★★★★★★★Outcome:', outcome_trait, '★★★★★★★'))
	Harmonisation(exposure_list, outcome_trait)
}




