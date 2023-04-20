#######################################################################
#             Biomarker and NAFLD and its clinical spectrum           #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS
# Last modified 13 Sep 2022


# This script conduct LDSC analysis, Pre-installed ldsc environment required, then source activate ldsc

#!/bin/bash
#@author: Dongze Chen


trait_list=( "Ala" "His" "Lac" "HDL" "TG" "Ferritin" "SerumIron" "TSAT" "VD" "IL6" "IL1RA" "NAFLD_Anstee" "NAFLD_FinnGen" "NAFLD_UKB" ) 
for trait in "${trait_list[@]}"
do
	./munge_sumstats.py --sumstats /mnt/data1/user/chendongze/Project/Biomarkers_NAFLD_Project/disease_${trait}.txt   --chunksize 500000 --out ./sumstats_data/Biomarkers_NAFLD_Project/disease_${trait} --merge-alleles w_hm3.snplist
done


exp_list=( "Ala" "His" "Lac" "HDL" "TG" "Ferritin" "SerumIron" "TSAT" "VD" "IL1RA" )
out_list=( "NAFLD_Anstee" "NAFLD_FinnGen" "NAFLD_UKB" )

for exp_trait in "${exp_list[@]}"
do
	for out_trait in "${out_list[@]}"
	do 
		./ldsc.py --rg ./sumstats_data/Biomarkers_NAFLD_Project/disease_${exp_trait}.sumstats.gz,./sumstats_data/Biomarkers_NAFLD_Project/disease_${out_trait}.sumstats.gz --ref-ld-chr eur_w_ld_chr/  --w-ld-chr eur_w_ld_chr/ --out ./GeneCorrResults/Biomarkers_NAFLD_Project/${exp_trait}_${out_trait}
	done
done


