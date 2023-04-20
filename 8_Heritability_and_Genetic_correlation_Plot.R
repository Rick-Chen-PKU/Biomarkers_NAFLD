#######################################################################
#             Biomarker and NAFLD and its clinical spectrum           #
#######################################################################
# > sessionInfo()
# R version 3.6.1 (2019-07-05)
# platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.4 LTS

# This script plot heritability estimation plot and genetic correlation plot


#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("D:/WorkDir/Rstudio/GWAS/Heritability_Genetic_correlation")
library(data.table)
library(ggplot2)
library(tidyverse)
library(corrplot)

if (!require("pacman")) install.packages("pacman")
pacman::p_load("data.table", "ggplot2", "tidyverse", "corrplot")



#---------------------------------------------------------------------#
#              heritability and genetic correlation                   #----
#---------------------------------------------------------------------#


heri <- readxl::read_excel("./Heritability.xlsx")
gc <- readxl::read_excel("./LDSC_summary.xlsx")

## Heritability ##
cairo_pdf('./Heritability.pdf', width = 5.5, height = 5)
heri %>% mutate(exposure = factor(exposure, levels = heri$exposure)) %>% mutate(h2 = round(h2, 2)) %>% 
  ggplot(mapping = aes(x = exposure, y = h2, fill = exposure)) + 
  geom_bar(stat = 'identity', 
           position = 'stack',   
           width = 0.8,
           color="black", linewidth = 0.3) +   
  geom_text(aes(label=h2), size = 3, position = position_stack(vjust = .5)) +
  geom_errorbar(aes(ymin=(h2),ymax=(h2+h2_se)), width=0.15, linewidth=0.3, position = "dodge", alpha=0.5)+
  labs(y = 'SNP heritability') +  
  xlab(NULL) +   
  theme_set(theme_bw()) +  
  theme_classic() +  
  scale_x_discrete(expand=c(0,1), position='bottom') +
  scale_y_continuous(expand = c(0, 0.005), breaks = c(0.11, 0.22, 0.33), limits = c(0, 0.33)) +  
  scale_fill_manual(values=c("#FABB6E", "#FC8002", "#ADDB88", "#369F2D", "#FAC7B3", "#AF7EC0", "#dad691",
                              "#b1c043", "#92C2DD", "#B0A875", "#63b2ee", "#dc96a1", "#8481BA", "#614099",
                              "#AF7EC0", "#B0A875", "#A42DEB", "#adb6b6", "#63b2ee")) +    
  theme(legend.text = element_text(size = 6)) + 
  theme(legend.position="none")+
  guides(fill=guide_legend(title.theme=element_text(size=6, face="bold"))) +
  theme(legend.key.size= unit(.2,"inches"),  
        legend.key.height= unit(.1,"inches"),
        legend.key.width= unit(.1,"inches") ) + 
  theme(text = element_text(size=10, face="bold"))+
  theme(axis.text.x = element_text(color = "black",size = 10, angle = 45, vjust = 0.6, face="bold"),
        axis.text.y = element_text(color = "black",size = 10, face="bold"))
dev.off()
 
## Genetic correlation  ## 
cairo_pdf('./Cor_Matrix.pdf', width = 7.42, height = 3)
col_set <- colorRampPalette(c("#77C034","white" ,"lightskyblue"), alpha = TRUE)
cor_Matrix <- gc %>% select(exposure, outcome, rg) %>% tidyr::spread(key = outcome, value = rg) %>% select(!exposure) %>%
  purrr::set_names("Anstee", "FinnGen", "UKB") %>% as.matrix() %>% t()
colnames(cor_Matrix) = unique(gc$exposure)
cor_Matrix %>% corrplot(method="square", type = "full", order = 'original', col = col_set(50), cl.pos = "r",
           tl.cex = 0.7, tl.col = "black", tl.srt = 20, tl.pos = "lt",
           addgrid.col = "grey70", outline ="grey60", cl.length=5,
           sig.level = .05, addCoef.col = "black", number.digits = 2, number.font = 1, number.cex = 0.8,
           insig = "n")
dev.off()

