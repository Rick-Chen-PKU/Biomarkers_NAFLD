#######################################################################
#             Biomarker and NAFLD and its clinical spectrum           #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS

# This script plot forestplots for clinical spectra


#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("D:/WorkDir/Rstudio/GWAS/Forestplot_New")
library(meta)
library(readxl)
library(tibble)
library(writexl)
library(forplo)
library(grid)
library(forestploter)

if (!require("pacman")) install.packages("pacman")
pacman::p_load("meta", "metafor", "openxlsx", "readxl", "writexl", "tidyverse", "dplyr", "tibble", "ggforestplot", "cowplot")

#---------------------------------------------------------------------#
#                   Load two-sample MR results                        #----
#---------------------------------------------------------------------#
dat1 <- data.frame(read_excel('./Biomarkers_NAFLD_Project_MR_Result_NAFL.xlsx'))
dat1$Study <- 'Anstee'
dat1 <- subset(dat1, select = c('Study', 'Exposure', 'Outcome', 'Method', 'beta', 'SE', "p_value", "OR", "OR_CI"))
dat1$p_value <- as.numeric(dat1$p_value)

dat2 <- data.frame(read_excel('./Biomarkers_NAFLD_Project_MR_Result_NASH.xlsx'))
dat2$Study <- 'FinnGen'
dat2 <- subset(dat2, select = c('Study', 'Exposure', 'Outcome', 'Method', 'beta', 'SE', "p_value", "OR", "OR_CI"))
dat2$p_value <- as.numeric(dat2$p_value)

dat3 <- data.frame(read_excel('./Biomarkers_NAFLD_Project_MR_Result_LF.xlsx'))
dat3$Study <- 'UKB'
dat3 <- subset(dat3, select = c('Study', 'Exposure', 'Outcome', 'Method', 'beta', 'SE', "p_value", "OR", "OR_CI"))
dat3$p_value <- as.numeric(dat3$p_value)

dat4 <- data.frame(read_excel('./Biomarkers_NAFLD_Project_MR_Result_Cirrhosis.xlsx'))
dat4$Study <- 'UKB'
dat4 <- subset(dat4, select = c('Study', 'Exposure', 'Outcome', 'Method', 'beta', 'SE', "p_value", "OR", "OR_CI"))
dat4$p_value <- as.numeric(dat4$p_value)

dat5 <- data.frame(read_excel('./Biomarkers_NAFLD_Project_MR_Result_HCC.xlsx'))
dat5$Study <- 'UKB'
dat5 <- subset(dat5, select = c('Study', 'Exposure', 'Outcome', 'Method', 'beta', 'SE', "p_value", "OR", "OR_CI"))
dat5$p_value <- as.numeric(dat5$p_value)


ukb_var <-  c('1_AST','1_ALT','1_GGT','1_TBIL','1_DBIL','1_TP','1_ALB', "2_Ile", "5_POA", "7_VEa","8_GDF15")  

#---------------------------------------------------------------------#
#                       Forestplot for NAFL                           #----
#---------------------------------------------------------------------#

library(extrafont)
# font_import()
loadfonts(device = "win")
windowsFonts("Arial" = windowsFont("Arial"))  

plot_data = dat1
plot_data <- plot_data[!plot_data$Exposure %in% ukb_var, ]
plot_data <- plot_data[plot_data$Method %in% c("IVW raw"), ]
# plot_data <- plot_data[plot_data$Method %in% c("IVW raw", "Weighted mode", "MR-RAPS"), ]
plot_data$OR <- as.numeric(plot_data$OR)
plot_data$LOR <- as.numeric(str_extract(plot_data$OR_CI, "(?<=\\().*(?=,)")) 
plot_data$UOR <- as.numeric(str_extract(plot_data$OR_CI, "(?<=\\,).*(?=\\))"))
plot_data$Group <- as.numeric(str_extract(plot_data$Exposure, "^.*(?=\\_)")) 
plot_data$Exposure <- str_extract(plot_data$Exposure, "(?<=\\_).*$")
plot_data <- subset(plot_data, select = c("Group", "Exposure", "Method", "OR", "LOR", "UOR", "p_value"))
# round(plot_data$p_value,4)

pdf(file="./Forestplot_Subprocess_NAFL.pdf", width=3, height=7, fonts="Arial", pointsize=6, onefile = FALSE)
forplo(plot_data[,c("OR", "LOR", "UOR")],
       sort=F, 
       left.align=TRUE, 
       favorlabs=c('Lower risk','Higher risk'),
       add.arrow.left=TRUE, 
       add.arrow.right=TRUE, 
       row.labels= c("Alanine", "Glutamine", "Homocysteine", "Histidine", "Leucine", "Phenylalanine", "Tyrosine", "Valine",
                     "Hemoglobin A1C", "Lactate",
                     "Apolipoprotein AI", "Apolipoprotein B", "HDL-cholesterol", "Lipoprotein(a)", "Triglycerides",
                     "Arachidonic acid", "Docosahexaenoic acid", "Docosapentaenoic acid", "Linoleic acid",
                     "Omega-3 fatty acid", "Omega-6 fatty acid", "Stearic acid", 
                     "Ferritin", "Serum iron", "TIBC", "Transferrin saturation",
                     "Folate", "Vitamin B12", "Vitamin C", "Vitamin D", 
                     "Adiponectin", "C1qTNF1", "C1qTNF5", "Complement C4", "C reaction protein", "Galectin-3", "Ghrelin",
                     "IL-27", "IL-18", "IL-1RA",  "IL-6",
                     "Leptin", "MMP7", "Periostin", "Resistin", "E-selectin", "Leukocyte telomere length"),   # plot_data$Exposure
       groups = plot_data$Group,   
       grouplabs=  c('Amino acid','Glucose', "Lipids", "Fatty acid", "Iron homeostasis", "Nutrients", "Cytokines", "Aging related trait"),
       shade.col = '#FDF0E0',  # shade.col='gray','#FDF0E0'
       shade.alpha = 0.9, 
       shade.every=1, 
       ci.edge=T,  
       ci.lwd = 0.6,  
       char=20, 
       size = 0.6,  
       col = "darkorange", 
       # right.bar = T, 
       title='NAFL outcome',
       margin.left = 12,
       margin.top = 0,
       margin.bottom = 2,
       margin.right = 10, 
       save = F,  
       save.path = "./",  #
       save.type = "pdf",  
       save.height = 7,   
       save.width = 12)   
dev.off() 



#---------------------------------------------------------------------#
#                       Forestplot for NASH                           #----
#---------------------------------------------------------------------#
plot_data = dat2
plot_data <- plot_data[!plot_data$Exposure %in% ukb_var, ]
plot_data <- plot_data[plot_data$Method %in% c("IVW raw"), ]
# plot_data <- plot_data[plot_data$Method %in% c("IVW raw", "Weighted mode", "MR-RAPS"), ]
plot_data$OR <- as.numeric(plot_data$OR)
plot_data$LOR <- as.numeric(str_extract(plot_data$OR_CI, "(?<=\\().*(?=,)")) 
plot_data$UOR <- as.numeric(str_extract(plot_data$OR_CI, "(?<=\\,).*(?=\\))"))
plot_data$Group <- as.numeric(str_extract(plot_data$Exposure, "^.*(?=\\_)")) 
plot_data$Exposure <- str_extract(plot_data$Exposure, "(?<=\\_).*$")
plot_data <- subset(plot_data, select = c("Group", "Exposure", "Method", "OR", "LOR", "UOR", "p_value"))
# round(plot_data$p_value,4)


pdf(file="./Forestplot_Subprocess_NASH.pdf", width=3, height=7, fonts="Arial", pointsize=6, onefile = FALSE)
forplo(plot_data[,c("OR", "LOR", "UOR")],
       sort=F, 
       left.align=TRUE, 
       favorlabs=c('Lower risk','Higher risk'),
       add.arrow.left=TRUE, 
       add.arrow.right=TRUE, 
       row.labels= c("Alanine", "Glutamine", "Homocysteine", "Histidine", "Leucine", "Phenylalanine", "Tyrosine", "Valine",
                     "Hemoglobin A1C", "Lactate",
                     "Apolipoprotein AI", "Apolipoprotein B", "HDL-cholesterol", "Lipoprotein(a)", "Triglycerides",
                     "Arachidonic acid", "Docosahexaenoic acid", "Docosapentaenoic acid", "Linoleic acid",
                     "Omega-3 fatty acid", "Omega-6 fatty acid", "Stearic acid", 
                     "Ferritin", "Serum iron", "TIBC", "Transferrin saturation",
                     "Folate", "Vitamin B12", "Vitamin C", "Vitamin D", 
                     "Adiponectin", "C1qTNF1", "C1qTNF5", "Complement C4", "C reaction protein", "Galectin-3", "Ghrelin",
                     "IL-27", "IL-18", "IL-1RA",  "IL-6",
                     "Leptin", "MMP7", "Periostin", "Resistin", "E-selectin", "Leukocyte telomere length"),   # plot_data$Exposure
       groups = plot_data$Group,   
       grouplabs=  c('Amino acid','Glucose', "Lipids", "Fatty acid", "Iron homeostasis", "Nutrients", "Cytokines", "Aging related trait"),
       shade.col = '#FDF0E0',  # shade.col='gray','#FDF0E0'
       shade.alpha = 0.9, 
       shade.every=1, 
       ci.edge=T,  
       ci.lwd = 0.6,  
       char=20, 
       size = 0.6,  
       col = "darkorange", 
       # right.bar = T, 
       title='NAFL outcome',
       margin.left = 12,
       margin.top = 0,
       margin.bottom = 2,
       margin.right = 10, 
       save = F,  
       save.path = "./",  #
       save.type = "pdf",  
       save.height = 7,   
       save.width = 12)   
dev.off() 


#---------------------------------------------------------------------#
#                      Forestplot for LF                              #----
#---------------------------------------------------------------------#
plot_data = dat3
plot_data <- plot_data[!plot_data$Exposure %in% ukb_var, ]
plot_data <- plot_data[plot_data$Method %in% c("IVW raw"), ]
# plot_data <- plot_data[plot_data$Method %in% c("IVW raw", "Weighted mode", "MR-RAPS"), ]
plot_data$OR <- as.numeric(plot_data$OR)
plot_data$LOR <- as.numeric(str_extract(plot_data$OR_CI, "(?<=\\().*(?=,)"))  
plot_data$UOR <- as.numeric(str_extract(plot_data$OR_CI, "(?<=\\,).*(?=\\))"))
plot_data$Group <- as.numeric(str_extract(plot_data$Exposure, "^.*(?=\\_)")) 
plot_data$Exposure <- str_extract(plot_data$Exposure, "(?<=\\_).*$")
plot_data <- subset(plot_data, select = c("Group", "Exposure", "Method", "OR", "LOR", "UOR", "p_value"))
# round(plot_data$p_value,4)

pdf(file="./Forestplot_Subprocess_LF.pdf", width=3, height=7, fonts="Arial", pointsize=6, onefile = FALSE)
forplo(plot_data[,c("OR", "LOR", "UOR")],
       sort=F, 
       left.align=TRUE, 
       favorlabs=c('Lower risk','Higher risk'),
       add.arrow.left=TRUE, 
       add.arrow.right=TRUE, 
       row.labels= c("Alanine", "Glutamine", "Homocysteine", "Histidine", "Leucine", "Phenylalanine", "Tyrosine", "Valine",
                     "Hemoglobin A1C", "Lactate",
                     "Apolipoprotein AI", "Apolipoprotein B", "HDL-cholesterol", "Lipoprotein(a)", "Triglycerides",
                     "Arachidonic acid", "Docosahexaenoic acid", "Docosapentaenoic acid", "Linoleic acid",
                     "Omega-3 fatty acid", "Omega-6 fatty acid", "Stearic acid", 
                     "Ferritin", "Serum iron", "TIBC", "Transferrin saturation",
                     "Folate", "Vitamin B12", "Vitamin C", "Vitamin D", 
                     "Adiponectin", "C1qTNF1", "C1qTNF5", "Complement C4", "C reaction protein", "Galectin-3", "Ghrelin",
                     "IL-27", "IL-18", "IL-1RA",  "IL-6",
                     "Leptin", "MMP7", "Periostin", "Resistin", "E-selectin", "Leukocyte telomere length"),   # plot_data$Exposure
       groups = plot_data$Group,   
       grouplabs=  c('Amino acid','Glucose', "Lipids", "Fatty acid", "Iron homeostasis", "Nutrients", "Cytokines", "Aging related trait"),
       shade.col = '#FDF0E0',  # shade.col='gray','#FDF0E0'
       shade.alpha = 0.9, 
       shade.every=1, 
       ci.edge=T,  
       ci.lwd = 0.6,  
       char=20, 
       size = 0.6,  
       col = "darkorange", 
       # right.bar = T, 
       title='NAFL outcome',
       margin.left = 12,
       margin.top = 0,
       margin.bottom = 2,
       margin.right = 10, 
       save = F,  
       save.path = "./",  #
       save.type = "pdf",  
       save.height = 7,   
       save.width = 12)   
dev.off() 


#---------------------------------------------------------------------#
#                    Forestplot for Cirrhosis                         #----
#---------------------------------------------------------------------#
plot_data = dat4
plot_data <- plot_data[!plot_data$Exposure %in% ukb_var, ]
plot_data <- plot_data[plot_data$Method %in% c("IVW raw"), ]
# plot_data <- plot_data[plot_data$Method %in% c("IVW raw", "Weighted mode", "MR-RAPS"), ]
plot_data$OR <- as.numeric(plot_data$OR)
plot_data$LOR <- as.numeric(str_extract(plot_data$OR_CI, "(?<=\\().*(?=,)")) 
plot_data$UOR <- as.numeric(str_extract(plot_data$OR_CI, "(?<=\\,).*(?=\\))"))
plot_data$Group <- as.numeric(str_extract(plot_data$Exposure, "^.*(?=\\_)")) 
plot_data$Exposure <- str_extract(plot_data$Exposure, "(?<=\\_).*$")
plot_data <- subset(plot_data, select = c("Group", "Exposure", "Method", "OR", "LOR", "UOR", "p_value"))
# round(plot_data$p_value,4)


pdf(file="./Forestplot_Subprocess_Cirrhosis.pdf", width=3, height=7, fonts="Arial", pointsize=6, onefile = FALSE)
forplo(plot_data[,c("OR", "LOR", "UOR")],
       sort=F, 
       left.align=TRUE, 
       favorlabs=c('Lower risk','Higher risk'),
       add.arrow.left=TRUE, 
       add.arrow.right=TRUE, 
       row.labels= c("Alanine", "Glutamine", "Homocysteine", "Histidine", "Leucine", "Phenylalanine", "Tyrosine", "Valine",
                     "Hemoglobin A1C", "Lactate",
                     "Apolipoprotein AI", "Apolipoprotein B", "HDL-cholesterol", "Lipoprotein(a)", "Triglycerides",
                     "Arachidonic acid", "Docosahexaenoic acid", "Docosapentaenoic acid", "Linoleic acid",
                     "Omega-3 fatty acid", "Omega-6 fatty acid", "Stearic acid", 
                     "Ferritin", "Serum iron", "TIBC", "Transferrin saturation",
                     "Folate", "Vitamin B12", "Vitamin C", "Vitamin D", 
                     "Adiponectin", "C1qTNF1", "C1qTNF5", "Complement C4", "C reaction protein", "Galectin-3", "Ghrelin",
                     "IL-27", "IL-18", "IL-1RA",  "IL-6",
                     "Leptin", "MMP7", "Periostin", "Resistin", "E-selectin", "Leukocyte telomere length"),   # plot_data$Exposure
       groups = plot_data$Group,   
       grouplabs=  c('Amino acid','Glucose', "Lipids", "Fatty acid", "Iron homeostasis", "Nutrients", "Cytokines", "Aging related trait"),
       shade.col = '#FDF0E0',  # shade.col='gray','#FDF0E0'
       shade.alpha = 0.9, 
       shade.every=1, 
       ci.edge=T,  
       ci.lwd = 0.6,  
       char=20, 
       size = 0.6,  
       col = "darkorange", 
       # right.bar = T, 
       title='NAFL outcome',
       margin.left = 12,
       margin.top = 0,
       margin.bottom = 2,
       margin.right = 10, 
       save = F,  
       save.path = "./",  #
       save.type = "pdf",  
       save.height = 7,   
       save.width = 12)   
dev.off() 


#---------------------------------------------------------------------#
#                         Forestplot for HCC                          #----
#---------------------------------------------------------------------#
plot_data = dat5
plot_data <- plot_data[!plot_data$Exposure %in% ukb_var, ]
plot_data <- plot_data[plot_data$Method %in% c("IVW raw"), ]
# plot_data <- plot_data[plot_data$Method %in% c("IVW raw", "Weighted mode", "MR-RAPS"), ]
plot_data$OR <- as.numeric(plot_data$OR)
plot_data$LOR <- as.numeric(str_extract(plot_data$OR_CI, "(?<=\\().*(?=,)"))  
plot_data$UOR <- as.numeric(str_extract(plot_data$OR_CI, "(?<=\\,).*(?=\\))"))
plot_data$Group <- as.numeric(str_extract(plot_data$Exposure, "^.*(?=\\_)")) 
plot_data$Exposure <- str_extract(plot_data$Exposure, "(?<=\\_).*$")
plot_data <- subset(plot_data, select = c("Group", "Exposure", "Method", "OR", "LOR", "UOR", "p_value"))
# round(plot_data$p_value,4)


pdf(file="./Forestplot_Subprocess_HCC.pdf", width=3, height=7, fonts="Arial", pointsize=6, onefile = FALSE)
forplo(plot_data[,c("OR", "LOR", "UOR")],
       sort=F, 
       left.align=TRUE, 
       favorlabs=c('Lower risk','Higher risk'),
       add.arrow.left=TRUE, 
       add.arrow.right=TRUE, 
       row.labels= c("Alanine", "Glutamine", "Homocysteine", "Histidine", "Leucine", "Phenylalanine", "Tyrosine", "Valine",
                     "Hemoglobin A1C", "Lactate",
                     "Apolipoprotein AI", "Apolipoprotein B", "HDL-cholesterol", "Lipoprotein(a)", "Triglycerides",
                     "Arachidonic acid", "Docosahexaenoic acid", "Docosapentaenoic acid", "Linoleic acid",
                     "Omega-3 fatty acid", "Omega-6 fatty acid", "Stearic acid", 
                     "Ferritin", "Serum iron", "TIBC", "Transferrin saturation",
                     "Folate", "Vitamin B12", "Vitamin C", "Vitamin D", 
                     "Adiponectin", "C1qTNF1", "C1qTNF5", "Complement C4", "C reaction protein", "Galectin-3", "Ghrelin",
                     "IL-27", "IL-18", "IL-1RA",  "IL-6",
                     "Leptin", "MMP7", "Periostin", "Resistin", "E-selectin", "Leukocyte telomere length"),   # plot_data$Exposure
       groups = plot_data$Group,   
       grouplabs=  c('Amino acid','Glucose', "Lipids", "Fatty acid", "Iron homeostasis", "Nutrients", "Cytokines", "Aging related trait"),
       shade.col = '#FDF0E0',  # shade.col='gray','#FDF0E0'
       shade.alpha = 0.9, 
       shade.every=1, 
       ci.edge=T,  
       ci.lwd = 0.6,  
       char=20, 
       size = 0.6,  
       col = "darkorange", 
       # right.bar = T, 
       title='NAFL outcome',
       margin.left = 12,
       margin.top = 0,
       margin.bottom = 2,
       margin.right = 10, 
       save = F,  
       save.path = "./",  #
       save.type = "pdf",  
       save.height = 7,   
       save.width = 12)   
dev.off() 


