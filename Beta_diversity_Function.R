###################################################################
###################################################################
###################################################################
###     Functions Beta-diversity -> Created by Adrià         ######
###  Linkedin: https://www.linkedin.com/in/adria-cruells/    ######
###  ORCID: https://orcid.org/0000-0002-1179-7997            ######
###  Github:    Working on it                                ######
###################################################################
###################################################################

#### THINKing in rarefy not sum/percentatge


###################################################
###############Chargue libraries###################
###################################################


pacman::p_load(vegan,####Analysis package
               caret,#### Analysis package
               phyloseq,#### Filetype package
               readxl, #### Entry table package
               openxlsx, #### Entry table package
               dplyr, #### Table manage package
               ggplot2, #### Plot package
               ggrepel, #### Plot package
               microbial ####normalitzation
) 


###################################################
###############Phyloseq entry######################
###################################################



## Automatic entry

  #Enter directory of a created phytloseq( *:RData)

## Manual entry

# Sample data

# Otu table

# Taxmat table


###################################################
###############BETA    PARAMETERS##################
###################################################
# 
# Physeq1 <- readRDS(file = "C:/Users/Adrià/Desktop/Trabajo/Microbioma R/Clean_codes/Gut.rds") #####Put the phyloseq object with sample data
# Matrix_type <- "wunifrac" ####### Could be bray, euclidean for low alfa diversity try Weighted Unifrac, could be better having a low rare species
# Normalitzation <- "relative" ###'"relative","TMM","vst","log2".
# Adjformula <- ("Post_coital")
# Grup_variable <- "Post_coital"####### Name the variable of interest Group
# Results_directory <- paste0(getwd(),"/Gràfics/Vaginal_samples/Beta/") ###### Can use getwd if you setted the wd before, put "/" at final
# Analysis_name <- "Model4"
###################################################
###############Beta FUNCTIONS######################
###################################################


Beta_function <- function(Physeq,Adjformula,Normalitzation,Matrix_type,Grup_variable){
  
  
  relab_genera <- microbial::normalize(Physeq, method = Normalitzation)
  

  ### Create tue phyloseq tree 
  random_tree = ape::rtree(ntaxa(relab_genera), rooted=TRUE, tip.label=taxa_names(relab_genera))
  
  relab_genera = merge_phyloseq(relab_genera,random_tree)
  
  ord_clr <- phyloseq::ordinate(relab_genera, "DPCoA")
  
  
  
  PC <-  merge(data.frame(ord_clr$li),data.frame(sample_data(relab_genera)),by = "row.names")
  
  
  a <- ggplot(PC,aes(Axis1,Post_coital,fill = Post_coital))+geom_boxplot()+
    scale_fill_manual(values=c("#00A087FF", "#3C5488FF","orange")) +
    xlab("")+ylab("")+
    theme_minimal()+
    theme(legend.position = "none",axis.text = element_text(size=5))
  b <- ggplot(PC,aes(Post_coital,Axis2,fill =Post_coital))+geom_boxplot()+
    theme_minimal()+
    xlab("")+ylab("")+
    scale_fill_manual(values=c("#00A087FF", "#3C5488FF","orange")) +
    theme(legend.position = "none",axis.text = element_text(size=5))

  c <- phyloseq::plot_ordination(relab_genera, ord_clr, type="samples", color="Pacient") + 
    geom_point(size = 6)+
    theme_minimal()+
    #coord_fixed((PC$Axis2 / PC$Axis1)) +
    scale_color_manual(values=c("#00A087FF", "#3C5488FF","orange")) +
    stat_ellipse(aes(group = Post_coital), linetype = 4)+
    theme(legend.position = "none")
  
  
  
  library(patchwork)
  
  design <- "
  1111#
  22223
  22223
  22223
  22223
"
  
  Beta_plot <- a + c + b + 
    plot_layout(design = design)
  
  ggsave(plot = Beta_plot,paste0(Results_directory,"/",Analysis_name,"_Beta.png"), width = 20, height = 20, units = "cm")
  
  
  
  
  clr_dist_matrix <- phyloseq::distance(relab_genera, method = Matrix_type) 
  
  dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(relab_genera)$Post_coital)

  
 
  
  boxplot(dispr, main = "", xlab = "")
  
  
  permutest(dispr) ### Hay diferencia estadistica

  #ADONIS test

  sink(file = paste0(Results_directory,"/",Analysis_name,"_Beta.txt"))
  print(vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(relab_genera)$Post_coital+phyloseq::sample_data(relab_genera)$Menopause+
                         phyloseq::sample_data(relab_genera)$BMI + phyloseq::sample_data(relab_genera)$Age))
  sink()
  
  return(Beta_plot)

}



###################################################
########RUN BETA   (For non source users)##########
###################################################

# Beta_function(Physeq,Matrix_type)

