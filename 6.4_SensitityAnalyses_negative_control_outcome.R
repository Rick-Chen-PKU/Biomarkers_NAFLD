#######################################################################
#             Biomarker and NAFLD and its clinical spectrum           #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022


# This script conduct main two-sample MR for negative control outcome
# including IVW, Simple median, Weighted median,Weighted mode, MR-Egger, DIVW, MR-RAPS, MR-PRESSO methods

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory
setwd("/home/dongzechen/GWAS/04-twosampleMR")  # your_working_directory

if (!require("pacman")){install.packages("pacman")}
pacman::p_load("TwoSampleMR",  "mr.raps", "MendelianRandomization", "MRPRESSO", "data.table", 
				"stringr", "readxl", "writexl", "stringr")

#---------------------------------------------------------------------#
#                         Run two-sample MR                           #----
#---------------------------------------------------------------------#
Project = 'Biomarkers_NAFLD_Project'
path = paste0("/mnt/data1/user/chendongze/Project/", Project, "/disease_")


NCO = function(exposure_trait, outcome_trait){

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


	# Search the IEU OpenGWAS database, fuction "extract_outcome_data" automatically finds the proxies.
	if (outcome_trait == 'Skin colour'){
		outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'ukb-b-19560') 
		outcome_dat$samplesize.outcome = 456692  # sample size information can be found in IEU OpenGWAS database
	} else if (outcome_trait == 'Ease of skin tanning'){
		outcome_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = 'ukb-b-533')
		outcome_dat$samplesize.outcome = 453065  # sample size information can be found in IEU OpenGWAS database
	} 


	dat = TwoSampleMR::harmonise_data(exp_dat, outcome_dat, action=2); nrow(dat)  # defalut action=2 would remove more SNPs, Remove SNPs for being palindromic with intermediate allele frequencies
	dat = dat[dat$ambiguous == 'FALSE', ]; nrow(dat)  # remove ambiguous SNPs
	dat = TwoSampleMR::steiger_filtering(dat) # steiger_filtering, must contain sample_size column
	dat = dat[dat$steiger_dir == 'TRUE', ]; nrow(dat) 
	MRInputObject <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse =  dat$se.outcome, 
	exposure = exposure_trait, outcome = outcome_trait, snps = dat$SNP, effect_allele = dat$effect_allele.exposure, other_allele = dat$other_allele.exposure, eaf = dat$eaf.exposure)
	# IVW
	IVWObject1 <- MendelianRandomization::mr_ivw(MRInputObject,model= "default",robust = FALSE, penalized = FALSE, correl = FALSE, weights ="simple", psi = 0,distribution = "normal",alpha = 0.05)
	# IVWObject1 
	vec1 = c(IVWObject1@Exposure, IVWObject1@Outcome, IVWObject1@SNPs, 'IVW raw', round(IVWObject1@Estimate,3), round(exp(IVWObject1@Estimate),3), 
			paste0("(", as.character(round(IVWObject1@CILower,3)),",", as.character(round(IVWObject1@CIUpper,3)), ")"),
			paste0("(", as.character(round(exp(IVWObject1@CILower),3)),",", as.character(round(exp(IVWObject1@CIUpper),3)), ")"), round(IVWObject1@StdError,3), IVWObject1@Pvalue,  round(IVWObject1@Heter.Stat[1],2), IVWObject1@Heter.Stat[2],
			round((IVWObject1@Heter.Stat[1] - IVWObject1@SNPs + 1 )/IVWObject1@Heter.Stat[1], 3) )

	# Mode-based estimation
	MBEObject <- mr_mbe(MRInputObject, weighting = "weighted", stderror = "delta", phi = 1, seed = 314159265, iterations = 10000, distribution = "normal", alpha = 0.05)
	# MBEObject
	vec4 = 	c(MBEObject@Exposure, MBEObject@Outcome, MBEObject@SNPs, 'Weighted mode', round(MBEObject@Estimate,3), round(exp(MBEObject@Estimate),3),
			paste0("(", as.character(round(MBEObject@CILower,3)),",", as.character(round(MBEObject@CIUpper,3)), ")"), 
			paste0("(", as.character(round(exp(MBEObject@CILower),3)),",", as.character(round(exp(MBEObject@CIUpper),3)), ")"), round(MBEObject@StdError,3), MBEObject@Pvalue, NA, NA, NA )


	# MR-RAPS
	res2 = mr.raps.mle(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, over.dispersion=TRUE, diagnostics=FALSE)
	vec7 = c(exposure_trait, outcome_trait, IVWObject1@SNPs, 'MR-RAPS', round(res2$beta.hat,3), round(exp(res2$beta.hat),3), 
				paste0('(', as.character(round(res2$beta.hat-1.96*res2$beta.se, 2)) , ',', as.character(round(res2$beta.hat+1.96*res2$beta.se, 2)), ')'), 
				paste0('(', as.character(round(exp(res2$beta.hat-1.96*res2$beta.se), 2)) , ',', as.character(round(exp(res2$beta.hat+1.96*res2$beta.se), 2)), ')'), 
				round(res2$beta.se,3), as.numeric(res2$beta.p.value), NA, NA, NA)


	out_df = rbind(vec1, vec4, vec7)
	out_df = data.frame(out_df)
	names(out_df) = c('Exposure', 'Outcome', 'N_snp', 'Method', 'beta', 'OR', 'beta_CI', 'OR_CI', 'SE', 'p_value', 'Q_Stat', 'Q_pvalue', 'I_square')

	# sensitivity analysis
	sen2 = mr_pleiotropy_test(dat)  # Horizontal pleiotropy test
	out_df$egger_intercept = sen2$egger_intercept
	out_df$egger_intercept_p_value = sen2$pval


	if (file.exists(paste0('./', Project, '_MR_Result_NCO.xlsx'))){
		old_dat = read_excel(paste0('./', Project, '_MR_Result_NCO.xlsx'))
		new_dat = rbind(old_dat,"",out_df)
		write_xlsx(new_dat, paste0('./', Project, '_MR_Result_NCO.xlsx')) 
	} else {
		write_xlsx(out_df, paste0('./', Project, '_MR_Result_NCO.xlsx')) 
	} 

}



exposure_list = list.files('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/NewVerIVs')
exposure_list = gsub(".csv", "", exposure_list)

outcome_list = c("Skin colour", "Ease of skin tanning")
for (exposure_trait in exposure_list){
	for (outcome_trait in outcome_list){
		NCO(exposure_trait, outcome_trait)
	}
}







