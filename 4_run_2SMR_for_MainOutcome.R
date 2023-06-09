#######################################################################
#             Biomarker and NAFLD and its clinical spectrum           #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS


# This script conduct main two-sample MR for main consortium outcome
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


TSMR <- function(exposure_list, outcome_trait){ 
	for (exposure_trait in exposure_list){
		dat = data.frame(fread(paste0('./Project/', Project, '/', exposure_trait, '_to_', outcome_trait, '_instruments.csv')))  
		if (sum(dat$mr_keep) >= 2) {
			MRInputObject <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse =  dat$se.outcome, 
			exposure = exposure_trait, outcome = outcome_trait, snps = dat$SNP, effect_allele = dat$effect_allele.exposure, other_allele = dat$other_allele.exposure, eaf = dat$eaf.exposure)
			# IVW
			IVWObject1 <- MendelianRandomization::mr_ivw(MRInputObject,model= "default",robust = FALSE, penalized = FALSE, correl = FALSE, weights ="simple", psi = 0,distribution = "normal",alpha = 0.05)
			# IVWObject1 
			vec1 = c(IVWObject1@Exposure, IVWObject1@Outcome, IVWObject1@SNPs, 'IVW raw', round(IVWObject1@Estimate,3), round(exp(IVWObject1@Estimate),3), 
					paste0("(", as.character(round(IVWObject1@CILower,3)),",", as.character(round(IVWObject1@CIUpper,3)), ")"),
					paste0("(", as.character(round(exp(IVWObject1@CILower),3)),",", as.character(round(exp(IVWObject1@CIUpper),3)), ")"), round(IVWObject1@StdError,3), IVWObject1@Pvalue,  round(IVWObject1@Heter.Stat[1],2), IVWObject1@Heter.Stat[2],
					round((IVWObject1@Heter.Stat[1] - IVWObject1@SNPs + 1 )/IVWObject1@Heter.Stat[1], 3) )
			# IVWObject2 <- MendelianRandomization::mr_ivw(MRInputObject,model= "default",robust = TRUE, penalized = TRUE, correl = FALSE, weights ="simple", psi = 0,distribution = "normal",alpha = 0.05)
			# IVWObject2 
			# vec2 = 	c(IVWObject2@Exposure, IVWObject2@Outcome, IVWObject2@SNPs, 'IVW robust', round(IVWObject2@Estimate,3), round(exp(IVWObject2@Estimate),3),
					# paste0("(", as.character(round(IVWObject2@CILower,3)),",", as.character(round(IVWObject2@CIUpper,3)), ")"),
					# paste0("(", as.character(round(exp(IVWObject2@CILower),3)),",", as.character(round(exp(IVWObject2@CIUpper),3)), ")"), round(IVWObject2@StdError,3), IVWObject2@Pvalue, NA, NA, NA )
			
			# median-based
			if (nrow(dat) > 2){
			MedianObject1 <-MendelianRandomization::mr_median(MRInputObject,weighting = "simple",distribution ="normal",alpha = 0.05,iterations = 10000,seed = 314159265) 
			# MedianObject1 
			vec2 = 	c(MedianObject1@Exposure, MedianObject1@Outcome, MedianObject1@SNPs, 'Simple median', round(MedianObject1@Estimate,3), round(exp(MedianObject1@Estimate),3),
					paste0("(", as.character(round(MedianObject1@CILower,3)),",", as.character(round(MedianObject1@CIUpper,3)), ")"),
					paste0("(", as.character(round(exp(MedianObject1@CILower),3)),",", as.character(round(exp(MedianObject1@CIUpper),3)), ")"), round(MedianObject1@StdError,3), MedianObject1@Pvalue, NA, NA, NA )
			MedianObject2 <-MendelianRandomization::mr_median(MRInputObject,weighting = "weighted",distribution ="normal",alpha = 0.05,iterations = 10000,seed = 314159265)  
			# MedianObject2
			vec3 = 	c(MedianObject2@Exposure, MedianObject2@Outcome, MedianObject2@SNPs, 'Weighted median', round(MedianObject2@Estimate,3), round(exp(MedianObject2@Estimate),3),
					paste0("(", as.character(round(MedianObject2@CILower,3)),",", as.character(round(MedianObject2@CIUpper,3)), ")"),
					paste0("(", as.character(round(exp(MedianObject2@CILower),3)),",", as.character(round(exp(MedianObject2@CIUpper),3)), ")"), round(MedianObject2@StdError,3), MedianObject2@Pvalue, NA, NA, NA )
			# MedianObject3 <-MendelianRandomization::mr_median(MRInputObject,weighting = "penalized",distribution ="normal",alpha = 0.05,iterations = 10000,seed = 314159265)  
			# MedianObject3
			} else {
				vec2 = c(exposure_trait, outcome_trait, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
				vec3 = c(exposure_trait, outcome_trait, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
			}
			
			# Mode-based estimation
			MBEObject <- mr_mbe(MRInputObject, weighting = "weighted", stderror = "delta", phi = 1, seed = 314159265, iterations = 10000, distribution = "normal", alpha = 0.05)
			# MBEObject
			vec4 = 	c(MBEObject@Exposure, MBEObject@Outcome, MBEObject@SNPs, 'Weighted mode', round(MBEObject@Estimate,3), round(exp(MBEObject@Estimate),3),
					paste0("(", as.character(round(MBEObject@CILower,3)),",", as.character(round(MBEObject@CIUpper,3)), ")"), 
					paste0("(", as.character(round(exp(MBEObject@CILower),3)),",", as.character(round(exp(MBEObject@CIUpper),3)), ")"), round(MBEObject@StdError,3), MBEObject@Pvalue, NA, NA, NA )
			
			# MR-Egger
			if (nrow(dat) > 2){
			EggerObject1 <-mr_egger(MRInputObject,robust = FALSE, penalized = FALSE,correl =FALSE,distribution = "normal",alpha = 0.05)
			# EggerObject1 
			vec5 = 	c(EggerObject1@Exposure, EggerObject1@Outcome, EggerObject1@SNPs, 'MR-Egger', round(EggerObject1@Estimate,3), round(exp(EggerObject1@Estimate),3),
					paste0("(", as.character(round(EggerObject1@CILower.Est,3)),",", as.character(round(EggerObject1@CIUpper.Est,3)), ")"), 
					paste0("(", as.character(round(exp(EggerObject1@CILower.Est),3)),",", as.character(round(exp(EggerObject1@CIUpper.Est),3)), ")"), round(EggerObject1@StdError.Est,3), EggerObject1@Pvalue.Est, NA, NA, NA )
			} else {
				vec5 = c(exposure_trait, outcome_trait, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
			}
			# EggerObject2 <-mr_egger(MRInputObject,robust = TRUE, penalized = TRUE,correl =FALSE,distribution = "normal",alpha = 0.05) 
			# EggerObject2 
			
			# Debiased inverse-variance weighted method
			# DIVWObject = mr_divw(MRInputObject, over.dispersion = TRUE, alpha = 0.05, diagnostics = FALSE)
			# DIVWObject
			# vec6 = 	c(DIVWObject@Exposure, DIVWObject@Outcome, DIVWObject@SNPs, 'DIVW', round(DIVWObject@Estimate,3), round(exp(DIVWObject@Estimate),3),
					# paste0("(", as.character(round(DIVWObject@CILower,3)),",", as.character(round(DIVWObject@CIUpper,3)), ")"),
					# paste0("(", as.character(round(exp(DIVWObject@CILower),3)),",", as.character(round(exp(DIVWObject@CIUpper),3)), ")"), round(DIVWObject@StdError,3), DIVWObject@Pvalue, NA, NA, NA )
			
			# MR-RAPS
			res2 = mr.raps.mle(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, over.dispersion=TRUE, diagnostics=FALSE)
			vec7 = c(exposure_trait, outcome_trait, IVWObject1@SNPs, 'MR-RAPS', round(res2$beta.hat,3), round(exp(res2$beta.hat),3), 
						paste0('(', as.character(round(res2$beta.hat-1.96*res2$beta.se, 2)) , ',', as.character(round(res2$beta.hat+1.96*res2$beta.se, 2)), ')'), 
						paste0('(', as.character(round(exp(res2$beta.hat-1.96*res2$beta.se), 2)) , ',', as.character(round(exp(res2$beta.hat+1.96*res2$beta.se), 2)), ')'), 
						round(res2$beta.se,3), as.numeric(res2$beta.p.value), NA, NA, NA)
		
			# MR-PRESSO
			if (nrow(dat) >= 4){
				res3 = mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
					OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 3000,  SignifThreshold = 0.05, seed=142857)
				vec81 = c(exposure_trait, outcome_trait, IVWObject1@SNPs,'MR-PRESSO:raw', round(res3$`Main MR results`$`Causal Estimate`[1],3), round(exp(res3$`Main MR results`$`Causal Estimate`[1]),2), 
					paste0('(', as.character(round(res3$`Main MR results`$`Causal Estimate`[1]-1.96*res3$`Main MR results`$`Sd`[1], 2)) , ',', as.character(round(res3$`Main MR results`$`Causal Estimate`[1]+1.96*res3$`Main MR results`$`Sd`[1], 2)), ')'), 
					paste0('(', as.character(round(exp(res3$`Main MR results`$`Causal Estimate`[1]-1.96*res3$`Main MR results`$`Sd`[1]), 2)) , ',', as.character(round(exp(res3$`Main MR results`$`Causal Estimate`[1]+1.96*res3$`Main MR results`$`Sd`[1]), 2)), ')'), 
					round(res3$`Main MR results`$`Sd`[1], 3), res3$`Main MR results`$`P-value`[1], NA, NA, NA)
				vec82 = c(exposure_trait, outcome_trait, IVWObject1@SNPs - length(res3$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),'MR-PRESSO:Outlier-corrected', 
					round(res3$`Main MR results`$`Causal Estimate`[2],3), round(exp(res3$`Main MR results`$`Causal Estimate`[2]),2), 
					paste0('(', as.character(round(res3$`Main MR results`$`Causal Estimate`[2]-1.96*res3$`Main MR results`$`Sd`[2], 2)) , ',', as.character(round(res3$`Main MR results`$`Causal Estimate`[2]+1.96*res3$`Main MR results`$`Sd`[2], 2)), ')'), 
					paste0('(', as.character(round(exp(res3$`Main MR results`$`Causal Estimate`[2]-1.96*res3$`Main MR results`$`Sd`[2]), 2)) , ',', as.character(round(exp(res3$`Main MR results`$`Causal Estimate`[2]+1.96*res3$`Main MR results`$`Sd`[2]), 2)), ')'), 
					round(res3$`Main MR results`$`Sd`[2], 3), res3$`Main MR results`$`P-value`[2], NA, NA, NA)
			} else {
				vec81 = c(exposure_trait, outcome_trait, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
				vec82 = c(exposure_trait, outcome_trait, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
			}


			out_df = rbind(vec1, vec2, vec3, vec4, vec5, vec7, vec81, vec82)
			out_df = data.frame(out_df)
			names(out_df) = c('Exposure', 'Outcome', 'N_snp', 'Method', 'beta', 'OR', 'beta_CI', 'OR_CI', 'SE', 'p_value', 'Q_Stat', 'Q_pvalue', 'I_square')
			
			# sensitivity analysis
			sen2 = mr_pleiotropy_test(dat)  # Horizontal pleiotropy test
			out_df$egger_intercept = sen2$egger_intercept
			out_df$egger_intercept_p_value = sen2$pval

			
			if (file.exists(paste0('./', Project, '_MR_Result_Anstee.xlsx'))){
				old_dat = read_excel(paste0('./', Project, '_MR_Result_Anstee.xlsx'))
				new_dat = rbind(old_dat,"",out_df)
				write_xlsx(new_dat, paste0('./', Project, '_MR_Result_Anstee.xlsx')) 
			} else {
				write_xlsx(out_df, paste0('./', Project, '_MR_Result_Anstee.xlsx')) 
			} 
			
			print(paste0('################################# Finish ', exposure_trait, ' to ', outcome_trait, ' MR analysis##################################'))
			print(paste0('################################# Finish ', exposure_trait, ' to ', outcome_trait, ' MR analysis##################################'))
		} else {
			print(paste0('################################# 暴露性状', exposure_trait, ' mr_keep个数不够 ##################################'))
		}
		
	}
}


# exposure_list = list.files('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/NewVerIVs')
# exposure_list = gsub(".csv", "", exposure_list)

exposure_list = c("4_TG2021")
outcome_list = c('NAFLD_Anstee', 'NAFLD_FinnGen', 'NAFLD_UKB')

for (outcome_trait in outcome_list){
	TSMR(exposure_list, outcome_trait)
}




