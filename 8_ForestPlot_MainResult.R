#######################################################################
#             Biomarker and NAFLD and its clinical spectrum           #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS

# This script plot the forestplots for NAFLD including FinnGen, Anstee, UKB and meta results

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

if (!require("pacman")) install.packages("pacman")
pacman::p_load("meta", "metafor", "openxlsx", "readxl", "writexl", "tidyverse", "dplyr", "tibble", "ggforestplot", "cowplot")

#---------------------------------------------------------------------#
#                   Load two-sample MR results                        #----
#---------------------------------------------------------------------#
dat1 <- data.frame(read_excel('./Biomarkers_NAFLD_Project_MR_Result_Anstee.xlsx'))
dat1$Study <- 'Anstee'
dat1 <- subset(dat1, select = c('Study', 'Exposure', 'Outcome', 'Method', 'beta', 'SE', "p_value", "OR", "OR_CI"))
dat1$p_value <- as.numeric(dat1$p_value)

dat2 <- data.frame(read_excel('./Biomarkers_NAFLD_Project_MR_Result_FinnGen.xlsx'))
dat2$Study <- 'FinnGen'
dat2 <- subset(dat2, select = c('Study', 'Exposure', 'Outcome', 'Method', 'beta', 'SE', "p_value", "OR", "OR_CI"))
dat2$p_value <- as.numeric(dat2$p_value)

dat3 <- data.frame(read_excel('./Biomarkers_NAFLD_Project_MR_Result_UKB.xlsx'))
dat3$Study <- 'UKB'
dat3 <- subset(dat3, select = c('Study', 'Exposure', 'Outcome', 'Method', 'beta', 'SE', "p_value", "OR", "OR_CI"))
dat3$p_value <- as.numeric(dat3$p_value)

data <- rbind(dat1, dat2, dat3)
data <- data.frame(data)
data <- na.omit(data)

ukb_var <-  c('1_AST','1_ALT','1_GGT','1_TBIL','1_DBIL','1_TP','1_ALB', "2_Ile", "5_POA", "7_VEa","8_GDF15") 
message(paste0('Total ', length(ukb_var), ' variables were derived from the UKB sample'))

data <- data[-which(data$Study == 'UKB' & data$Exposure %in% ukb_var), ]; rownames(data) <- NULL

test_exp <- 'TP'
dat1[which(dat1$Exposure==test_exp & dat1$Method=='IVW raw'), 'beta']
dat2[which(dat2$Exposure==test_exp & dat2$Method=='IVW raw'), 'beta']
dat3[which(dat3$Exposure==test_exp & dat3$Method=='IVW raw'), 'beta']

#---------------------------------------------------------------------#
#                           Run meta-analysis                         #----
#---------------------------------------------------------------------#
# Function to run a fixed effect meta-analysis 
meta_func <- function(method_varname, exp_varname){
  temp_data <- data[data$Exposure == exp_varname & data$Method==method_varname,]
  temp_data$beta = as.numeric(temp_data$beta)
  temp_data$SE = as.numeric(temp_data$SE)
  m <- metagen(beta, SE, data=temp_data, studlab=paste(Study), fixed = TRUE, random = FALSE, prediction=TRUE, sm="SMD") 
  #extract values from meta output
  TE.tibble <- as_tibble(m$TE.fixed)
  se.tibble <- as_tibble(m$seTE.fixed)
  p.tibble <- as_tibble(m$pval.fixed)
  #combine tibbles and change column names
  tibble <- cbind(TE.tibble, se.tibble, p.tibble)
  colnames(tibble) <- c("b", "se", "pval")
  # add columns for exposure, outcome and method
  tibble$exposure <- exp_varname
  tibble$outcome <- paste0(temp_data$Study, collapse="-")
  tibble$method <- method_varname
  return (data.frame(tibble))
}

# unique(data$Method)
# method_list = c("IVW raw", "Simple median", "Weighted median", "Weighted mode", "MR-Egger", "MR-RAPS", "MR-PRESSO:raw")
# exp_list = unique(data$Exposure)
# out_list = unique(data$Outcome)

res <- data.frame()
var_pair <- unique(data[,c('Method', 'Exposure')])

for (i in 1:nrow(var_pair )){
  method_varname = var_pair[i, 'Method']
  exp_varname = var_pair[i, 'Exposure']
  tp <- meta_func(method_varname, exp_varname)
  res <- rbind(res, tp)
}
write_xlsx(res, './meta_analysis_res.xlsx')


#---------------------------------------------------------------------#
#                     Forestplot for Anstee                           #----
#---------------------------------------------------------------------#
library(extrafont)
# font_import()
loadfonts(device = "win")
windowsFonts("Arial" = windowsFont("Arial")) 

plot_data = dat1
ukb_var <-  c('1_AST','1_ALT','1_GGT','1_TBIL','1_DBIL','1_TP','1_ALB', "2_Ile", "5_POA", "7_VEa","8_GDF15")  
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

# tiff(file="./Figure/Forestplot_Univariate.tif", width=6, height=20, units = 'in', res=600, pointsize=8)
pdf(file="./Forestplot_Anstee.pdf", width=3, height=7, fonts="Arial", pointsize=6, onefile = FALSE)
forplo(plot_data[,c("OR", "LOR", "UOR")],
       # font = 'Arial',  
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
       title='Forestplot from Anstee outcome',
       margin.left = 12,
       margin.top = 0,
       margin.bottom = 2,
       margin.right = 10, 
       save = F,  
       save.path = "./",  
       # save.name = "Forestplot of multivariable Logistic regression model",
       save.type = "pdf",  
       save.height = 7,  
       save.width = 12)   
dev.off() 


#---------------------------------------------------------------------#
#                     Forestplot for FinnGen                          #----
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


# tiff(file="./Figure/Forestplot_Univariate.tif", width=6, height=20, units = 'in', res=600, pointsize=8)
pdf(file="./Forestplot_FinnGen.pdf", width=3, height=7, fonts="Arial", pointsize=6, onefile = FALSE)
forplo(plot_data[,c("OR", "LOR", "UOR")],
       # font = 'Arial',  
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
       title='Forestplot from Anstee outcome',
       margin.left = 12,
       margin.top = 0,
       margin.bottom = 2,
       margin.right = 10, 
       save = F,  
       save.path = "./",  
       # save.name = "Forestplot of multivariable Logistic regression model",
       save.type = "pdf",  
       save.height = 7,  
       save.width = 12)   
dev.off() 


#---------------------------------------------------------------------#
#                     Forestplot for UKB                              #----
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


# tiff(file="./Figure/Forestplot_Univariate.tif", width=6, height=20, units = 'in', res=600, pointsize=8)
pdf(file="./Forestplot_UKB.pdf", width=3, height=7, fonts="Arial", pointsize=6, onefile = FALSE)
forplo(plot_data[,c("OR", "LOR", "UOR")],
       # font = 'Arial',  
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
       title='Forestplot from Anstee outcome',
       margin.left = 12,
       margin.top = 0,
       margin.bottom = 2,
       margin.right = 10, 
       save = F,  
       save.path = "./",  
       # save.name = "Forestplot of multivariable Logistic regression model",
       save.type = "pdf",  
       save.height = 7,  
       save.width = 12)   
dev.off() 


#---------------------------------------------------------------------#
#                  Forestplot for Meta-analysis                       #----
#---------------------------------------------------------------------#
plot_data = data.frame(read_excel('./meta_analysis_res.xlsx'))
names(plot_data) <- c("b", "se", "p_value","Exposure", "Outcome", "Method")
plot_data <- plot_data[!plot_data$Exposure %in% ukb_var, ]
plot_data <- plot_data[plot_data$Method %in% c("IVW raw"), ]
# plot_data <- plot_data[plot_data$Method %in% c("IVW raw", "Weighted mode", "MR-RAPS"), ]
plot_data$OR <- as.numeric(exp(plot_data$b))
plot_data$LOR <- as.numeric(exp(plot_data$b-1.96*plot_data$se)) 
plot_data$UOR <- as.numeric(exp(plot_data$b+1.96*plot_data$se))
plot_data$Group <- as.numeric(str_extract(plot_data$Exposure, "^.*(?=\\_)")) 
plot_data$Exposure <- str_extract(plot_data$Exposure, "(?<=\\_).*$")
plot_data <- subset(plot_data, select = c("Group", "Exposure", "Method", "OR", "LOR", "UOR", "p_value"))
# round(plot_data$p_value,4)


# tiff(file="./Figure/Forestplot_Univariate.tif", width=6, height=20, units = 'in', res=600, pointsize=8)
pdf(file="./Forestplot_Meta.pdf", width=3, height=7, fonts="Arial", pointsize=6, onefile = FALSE)
forplo(plot_data[,c("OR", "LOR", "UOR")],
       # font = 'Arial',  
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
       title='Forestplot from Anstee outcome',
       margin.left = 12,
       margin.top = 0,
       margin.bottom = 2,
       margin.right = 10, 
       save = F,  
       save.path = "./",  
       # save.name = "Forestplot of multivariable Logistic regression model",
       save.type = "pdf",  
       save.height = 7,  
       save.width = 12)   
dev.off() 



