#######################################################################
#             Biomarker and NAFLD and its clinical spectrum           #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS

# This script plot forestplot for MVMR results


#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("D:/WorkDir/Rstudio/GWAS/Forestplot_MVMR")
library(readxl)
library(tibble)
library(writexl)
library(forplo)
library(tidyverse)
library(grid)
library(forestploter)

if (!require("pacman")) install.packages("pacman")
pacman::p_load("meta", "metafor", "openxlsx", "readxl", "writexl", "tidyverse", "dplyr", "tibble", "ggforestplot", "cowplot")

#---------------------------------------------------------------------#
#                   Load two-sample MR results                        #----
#---------------------------------------------------------------------#

mvmrda <- read_excel('./plot_data_manual.xlsx') %>% 
  mutate(or = exp(beta), lor = exp(beta - 1.96*se), uor = exp(beta + 1.96*se)) %>% data.frame()



#---------------------------------------------------------------------#
#                       Forestplot for MVMR                           #----
#---------------------------------------------------------------------#

library(extrafont)
# font_import()
loadfonts(device = "win")
windowsFonts("Arial" = windowsFont("Arial")) 
pdf(file="./Forestplot_MVMR.pdf", width=5, height=5, fonts="Arial", pointsize=6, onefile = FALSE)
forplo(mvmrda[,c("or", "lor", "uor")],
       sort=F, 
       left.align=TRUE,
       favorlabs=c('Lower risk','Higher risk'),
       add.arrow.left=TRUE, 
       add.arrow.right=TRUE, 
       row.labels= mvmrda$adjust,   # plot_data$Exposure
       groups = c(rep(1,5),rep(2,5),rep(3,5),rep(4,5),rep(5,5),rep(6,5),rep(7,5)),  
       grouplabs=  c('MR study to evaluate the causal association of Ala on NAFLD',
                     'MR study to evaluate the causal association of HDL on NAFLD', 
                     "MR study to evaluate the causal association of TG on NAFLD", 
                     "MR study to evaluate the causal association of Ferritin on NAFLD", 
                     "MR study to evaluate the causal association of SerumIron on NAFLD",
                     "MR study to evaluate the causal association of TSAT on NAFLD", 
                     "MR study to evaluate the causal association of VD on NAFLD"),
       shade.col = '#FDF0E0',  # shade.col='gray','#FDF0E0'
       shade.alpha = 0.9, 
       shade.every=1,  
       ci.edge=T,  
       ci.lwd = 0.6,  
       char=20,  
       size = 0.6,  
       col = "darkorange",  
       # right.bar = T, 
       title='NAFLD_FinnGen as outcome',
       margin.left = 28,
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






