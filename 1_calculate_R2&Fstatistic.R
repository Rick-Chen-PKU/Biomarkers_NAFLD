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
# This script estimates R2 and F-statistics for each of the exposures


#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls())   #  Remove any existing objects in R 

# Set working directory 
setwd('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/NewVerIVs/')   # your_working_directory


#---------------------------------------------------------------------#
#                            Reading Exposure                         #----
#---------------------------------------------------------------------#
# IVs are pre-extracted directly from the source literature and saved in csv format.
library(data.table)
library(tidyverse)
library(writexl)

if (file.exists('exp_data.xlsx')){
	system('rm exp_data.xlsx')
}


file_list = list.files()
i = 1
exp_dat = data.frame()  # for combine exposures
while (i <= length(file_list)){
	message(i)
	assign(paste0('exp_dat_', i), fread(file_list[i]))
	if (i >= 2){
		if (nrow(exp_dat)==0){
			exp_dat_1 = subset(exp_dat_1, select=c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'N', 'Trait'))
			exp_dat_2 = subset(exp_dat_2, select=c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'N', 'Trait'))
			exp_dat = rbind(exp_dat_1, exp_dat_2)
		} else {
			tmp_dat = get(paste0('exp_dat_', i))
			tmp_dat = subset(tmp_dat, select=c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'N', 'Trait'))
			exp_dat = rbind(exp_dat, tmp_dat)
		}
	}
	i = i + 1
}


#---------------------------------------------------------------------#
#                          R2 and F-statistic                         #----
#---------------------------------------------------------------------#

# Calculate R2 and F statistics for each exposure dataset
# method 1
# exp_dat$MAF <- ifelse(exp_dat$EAF>=0.5, 1-exp_dat$EAF, exp_dat$EAF)
# exp_dat$r2 <- 2 * (abs(exp_dat$BETA)) ** 2 * exp_dat$MAF * (1 - exp_dat$MAF) 
# exp_dat$F <- exp_dat$r2 * (exp_dat$N - 2) / (1 - exp_dat$r2)

# method 2
# exp_dat$r2 <- (2 * (exp_dat$BETA^2) * exp_dat$EAF * (1 - exp_dat$EAF)) /
  # (2 * (exp_dat$BETA^2) * exp_dat$EAF * (1 - exp_dat$EAF) + 2 * exp_dat$N * exp_dat$EAF * 
	# (1 - exp_dat$EAF) * exp_dat$SE^2)
# exp_dat$F <- exp_dat$r2 * (exp_dat$N - 2) / (1 - exp_dat$r2)

# method 3
exp_dat$F <- exp_dat$BETA^2 / exp_dat$SE^2
exp_dat$r2 <- exp_dat$F/(exp_dat$N-2+exp_dat$F)

# Calculate total R2 and F-statistic for each exposure dataset 
df1 = aggregate(exp_dat[,'r2'], by=list(Trait = exp_dat$Trait), FUN=sum)
df1$r2 = df1$r2 * 100
df2 = aggregate(exp_dat[,'F'], by=list(Trait = exp_dat$Trait), FUN=min)
names(df2) = c('Trait', 'F_min')
df3 = aggregate(exp_dat[,'F'], by=list(Trait = exp_dat$Trait), FUN=max)
names(df3) = c('Trait', 'F_max')

data_list <- list(df1, df2, df3)
df <- data_list %>% reduce(inner_join, by = "Trait")

# save
write_xlsx(exp_dat, "/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/exp_data.xlsx")
write_xlsx(df, "/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/Biomarkers_NAFLD_Project/exp_data_aggr.xlsx")




