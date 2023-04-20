# Biomarkers_NAFLD
This code relates to the project titled "Association of 60 circulating biomarkers with non-alcoholic fatty liver disease and its clinical spectrum: an integrative Mendelian randomization study of European ancestry".

# Step 0: Data preparation and pre-processing
Prepare the required GWAS summary data or extract the IV directly from the source literature

0_extract_IVs.R

# Step 1: Calculate R2 and F-statistic
First, we needed to calculate R2 and F-statistic for 60 biomarkers instrumental variables. We downloaded complete summary statistics for online resources.

1_calculate_R2&Fstatistic.R

# Step 2: Find LD proxies for missing SNPs
Next, we needed to find linkage disequilibrium (LD) proxies for genetic variants missing in the outcome dataset. Here, we did this using the LDlink.

2_find_LD_proxies_for_MainConsortium.R


# Step 3: Harmonisation
This step conduct Harmonisation for and the IVs with different outcome, produce instrumental variables Output in two-sample MR format.

3_harmonisation.R


# Step 4: Run two-sample MR analysis
Once our exposure and outcome datasets were ready, we ran the two-sample MR analysis.

4_run_2SMR_for_MainOutcome.R  


# Step 5: Run meta-analysis
We combined two-sample MR results obtained in the previous step using a fixed-effect meta-analysis. This was done using the meta package.

5_Meta_analysis.R

# Step 6: Sensitity Analyses
This step focuses on performing sensitivity analysis, including "risk_factor_outcome" and "negative_control_outcome".

6.4_SensitityAnalyses_negative_control_outcome.R
6.5_SensitityAnalyses_risk_factor_outcome.R

# Step 7: Secondary Analyses
This step focuses on performing secondary analysis, including "Subprocess", "LDSC" and "MVMR".

6.1_find_LD_proxies_for_Subprocess.R
6.2_harmonisation_Subprocess.R
6.3_run_2SMR_for_Subprocess.R
7_SecondaryAnalyses_LDSC.sh 
7_SecondaryAnalyses_MVMR_with_function.R

# Step 8: Create visualisations
This step create graphics diagnostics visualisation for main outcome two-sample MR including scatter plot, funnel plot, forest plot and leave-one-out plot. We created forest plots for discovery, replication, meta analyses, MVMR and subprocess using the forplo packages. We also created heritability estimates and genetic association correlation coefficient plots.

8_Create_graphics_diagnostics_visualisations.R 
8_ForestPlot_MVMR.R
8_ForestPlot_MainResult.R
8_ForestPlot_Subprocess.R
8_Heritability_and_Genetic_correlation_Plot.R

