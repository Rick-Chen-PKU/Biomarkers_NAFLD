#######################################################################
#             Biomarker and NAFLD and its clinical spectrum           #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022

#---------------------------------------------------------------------#
# 						Usage Description							  #
#---------------------------------------------------------------------#
# This script finds LD-proxies for subprocess outcomes, which can then be used in two-sample MR analyses
# Output in two-sample MR format 

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory
setwd("/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project")  # your_working_directory

if (!require("pacman")){install.packages("pacman")}
pacman::p_load("TwoSampleMR", "tidyverse", "dplyr", "data.table",  "LDlinkR", "readxl", "writexl", "stringr")

#---------------------------------------------------------------------#
#                 Read exposure and outcome datasets                  #----
#---------------------------------------------------------------------#

# Read exposure
exp_dat <- read_excel("exp_data.xlsx", col_names = T)

# Format exposure
if ('EAF' %in% colnames(exp_dat)){
		exp_dat = format_data(exp_dat,
			type = "exposure",
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
		exp_dat = format_data(exp_dat,
			type = "exposure",
			phenotype_col = "Trait",
			snp_col = "SNP",
			beta_col = "BETA",
			se_col = "SE",
			effect_allele_col = "A1",
			other_allele_col = "A2",
			pval_col = "P",
			samplesize_col = "N")
}

# Function to read outcome
out_func_practical <- function(filepath, var)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_var <- fread(paste0(filepath, 'disease_', var, '.txt'), na.strings="")
  outcome_var$BETA <- log(as.numeric(outcome_var$OR))
  outcome_var$Trait = var
  return(outcome_var)
}

filepath = '/mnt/data1/user/chendongze/Project/Biomarkers_NAFLD_Project/'  # your_summary_data_directory
out_nafl <- out_func_practical(filepath, "NAFL_UKBB")
out_nash <- out_func_practical(filepath, "NASH_FinnGen")
out_lf <- out_func_practical(filepath, "LF_FinnGen")
out_cirrhosis <- out_func_practical(filepath, "Cirrhosis_FinnGen")
out_hcc <- out_func_practical(filepath, "HCC_FinnGen")
out_cirrhosis_emdin <- out_func_practical(filepath, "Cirrhosis_Emdin")

#---------------------------------------------------------------------#
#              Identify SNPs that need proxies                        #----
#---------------------------------------------------------------------#

# Function to find list of snps in exposure dataset that are missing from the outcome dataset
find_missing_SNP <- function(out_dat) {
    snps_need_proxy <- setdiff(exp_dat$SNP, intersect(exp_dat$SNP, out_dat$SNP))
	return (snps_need_proxy)
}

# trait_list = unique(exp_dat$exposure)

mis_nafl <- find_missing_SNP(out_nafl) 
mis_nash <- find_missing_SNP(out_nash) 
mis_lf <- find_missing_SNP(out_lf) 
mis_cirrhosis <- find_missing_SNP(out_cirrhosis) 
mis_hcc <- find_missing_SNP(out_hcc) 
mis_cirrhosis_emdin <- find_missing_SNP(out_cirrhosis_emdin) 

# length(mis_nafl) # check how many snps are in vector
# mis_nafl[1] #see snp missing n#1
# mis_nafl[2] #see snp missing n#2
# mis_nafl[3] #see snp missing n#3

#---------------------------------------------------------------------#
#                    Find proxies for these snps                      #----
#---------------------------------------------------------------------#
# Function to find LD proxy using LDLINK
# We need to make sure the identified proxy SNPs are available in outcome dataset before continuing 
# If SNPs aren't available, you need to find the next best proxy for the missing SNP
# We exclude SNPs for which no suitable proxy can be found

my_token = "bf61a241c044"  # Need to request for LDlink VIP request (https://ldlink.nci.nih.gov/?tab=apiaccess)
find_LD_proxy <- function(snps_need_proxy, out_dat) {
		proxy_list = list()
		snp = out_dat$SNP
		for (i in 1:length(snps_need_proxy)){
			print(paste0('---### Total ', length(snps_need_proxy), ' SNPs, rank ', i, ' SNP ### ---'))
			system(paste0("curl -o result.txt -k -X GET 'https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=", snps_need_proxy[i],"&pop=EUR&r2_d=r2&window=20000&genome_build=grch37&token=bf61a241c044'"))
			# temp <- LDproxy(snps_need_proxy[i], "EUR", "r2", token = my_token, file = F, genome_build = "grch37")  # encounter Bad Gateway (HTTP 502)
			temp = data.frame(fread('./result.txt', quote=","))
			if (!TRUE %in% grepl('error', temp[,1])){
				proxy = as.character(temp[temp$R2 >= 0.8,'RS_Number'])
				temp1 <- intersect(proxy, snp)
				if (length(temp1) != 0){
					proxy_list[snps_need_proxy[i]] =  temp1[1]
				} else {
					proxy_list[snps_need_proxy[i]] =  '<<<Nope>>>'
				}
			} else {
				proxy_list[snps_need_proxy[i]] =  '<<<Nope>>>'
			}
			system("rm result.txt")
		}
		return (proxy_list)
		# we could change number of proxies we want to find
}


# These proxies imformation would be used in next section !!!
proxy_snp_nafl  <- find_LD_proxy(mis_nafl, out_nafl)
proxy_snp_1 = unlist(proxy_snp_nafl)
names(proxy_snp_1) = names(proxy_snp_nafl)

proxy_snp_nash  <- find_LD_proxy(mis_nash, out_nash)  
proxy_snp_2 = unlist(proxy_snp_nash)
names(proxy_snp_2) = names(proxy_snp_nash)

proxy_snp_lf  <- find_LD_proxy(mis_lf, out_lf)  
proxy_snp_3 = unlist(proxy_snp_lf)
names(proxy_snp_3) = names(proxy_snp_lf)

proxy_snp_cirrhosis  <- find_LD_proxy(mis_cirrhosis, out_cirrhosis)  
proxy_snp_4 = unlist(proxy_snp_cirrhosis)
names(proxy_snp_4) = names(proxy_snp_cirrhosis)

proxy_snp_hcc  <- find_LD_proxy(mis_hcc, out_hcc)
proxy_snp_5 = unlist(proxy_snp_hcc)
names(proxy_snp_5) = names(proxy_snp_hcc)

# save 
proxy_snp_1_df = data.frame(need_proxy=names(proxy_snp_1), proxy_snp=unname(proxy_snp_1))
proxy_snp_1_df = proxy_snp_1_df[proxy_snp_1_df$proxy_snp != '<<<Nope>>>', ]
writexl::write_xlsx(proxy_snp_1_df, './proxy_snp_nafl.xlsx')

proxy_snp_2_df = data.frame(need_proxy=names(proxy_snp_2), proxy_snp=unname(proxy_snp_2))
proxy_snp_2_df = proxy_snp_2_df[proxy_snp_2_df$proxy_snp != '<<<Nope>>>', ]
writexl::write_xlsx(proxy_snp_2_df, './proxy_snp_nash.xlsx')

proxy_snp_3_df = data.frame(need_proxy=names(proxy_snp_3), proxy_snp=unname(proxy_snp_3))
proxy_snp_3_df = proxy_snp_3_df[proxy_snp_3_df$proxy_snp != '<<<Nope>>>', ]
writexl::write_xlsx(proxy_snp_3_df, './proxy_snp_lf.xlsx')

proxy_snp_4_df = data.frame(need_proxy=names(proxy_snp_4), proxy_snp=unname(proxy_snp_4))
proxy_snp_4_df = proxy_snp_4_df[proxy_snp_4_df$proxy_snp != '<<<Nope>>>', ]
writexl::write_xlsx(proxy_snp_4_df, './proxy_snp_cirrhosis.xlsx')

proxy_snp_5_df = data.frame(need_proxy=names(proxy_snp_5), proxy_snp=unname(proxy_snp_5))
proxy_snp_5_df = proxy_snp_5_df[proxy_snp_5_df$proxy_snp != '<<<Nope>>>', ]
writexl::write_xlsx(proxy_snp_5_df, './proxy_snp_hcc.xlsx')

