###################################################################
###################################################################
###################################################################
###     Functions Main file exam -> Created by Adrià         ######
###  Linkedin: https://www.linkedin.com/in/adria-cruells/    ######
###  ORCID: https://orcid.org/0000-0002-1179-7997            ######
###  Github:    Working on it                                ######
###################################################################
###################################################################

rm(list=ls())

######## Enter the phyloseq object, set working directory and chargue libraries



setwd("C:/Users/Adrià/Desktop/Trabajo/Microbioma R") ###### put the directory were you work
Vaginal <- readRDS(file = "C:/Users/Adrià/Desktop/Trabajo/Microbioma R/Clean_codes/Vaginal.rds")  #Enter directory of a created phytloseq( *:RData)
Gut <- readRDS(file = "C:/Users/Adrià/Desktop/Trabajo/Microbioma R/Clean_codes/Gut.rds")  #Enter directory of a created phytloseq( *:RData)    subset_taxa(!Genus == "Homo")%>% for vaginal and fecal microbiome


##########
####Filter
##########
library(dplyr)
library(phyloseq)
library(microViz)

Vaginal <- Vaginal %>%
  subset_taxa(!grepl("unk",Genus))%>%
  subset_taxa(!grepl("Ralstonia",Genus))%>%
  subset_taxa(!grepl("unk",Species))

Gut <- Gut %>%
  subset_taxa(!grepl("unk",Genus))%>%
  subset_taxa(!grepl("unk",Species))



##########
####Create the variable for Controls // ITU_post coital // ITU no post-coital
##########
Gut <- Gut %>%
ps_mutate(
  Post_coital = case_when(Pacient == "No" ~ "Control",
                          Pacient == "Yes" &  Postcoital.UTI == "Yes"~ "UTI-postcoital",
                          Pacient == "Yes" &  Postcoital.UTI == "No"  ~ "UTI-no-postcoital")
)

Vaginal <- Vaginal %>%
  ps_mutate(
    Post_coital = case_when(Pacient == "No" ~ "Control",
                            Pacient == "Yes" &  Postcoital.UTI == "Yes"~ "UTI-postcoital",
                            Pacient == "Yes" &  Postcoital.UTI == "No"  ~ "UTI-no-postcoital")
  )
##########
#### Filter only the pattients to check the effect of pelvic score
##########



Vaginal <- Vaginal %>%
  ps_filter(Pacient == "Yes")

Gut <- Gut %>%
  ps_filter(Pacient == "Yes")



##########################
####Patients data analysis
##########################


#### obtain the data from phyloseq and select the interest variables ( in this case (Patient//Age(years)//BMI//Menopause//UTI/year))


# Patients_data <- data.frame(sample_data(Vaginal)) %>%
#   dplyr::select(Pacient2,BMI,Age,Menopause,UTI.year) 
# 
# Patients_data$UTI.year <- as.numeric(Patients_data$UTI.year)
# Patients_data$UTI.year[Patients_data$Pacient2 == "Control"] <- 0
# # Factor the basic variables that
# # we're interested in
# 
# ####Rename all the variables to get the correct table 
# 
# Patients_data$Menopause <- 
#   factor(Patients_data$Menopause, levels=c("Yes",
#                                            "No"),
#          labels=c("Yes", 
#                   "No"))
# 
# label(Patients_data$BMI)       <- "Body mass index (BMI)"
# label(Patients_data$Menopause)       <- "Menopause"
# label(Patients_data$Age)     <- "Age (Years)"
# label(Patients_data$UTI.year) <- "UTI / Year"
# 
# 
# library(table1)
# 
# 
# table1(~  BMI + Menopause + Age + UTI.year | Pacient2, data=Patients_data)

################
##Analysis parameters
################

Adjformula <- ("Score_pelvic+Age+BMI+Menopause") #####Model formula (Ex model 3)
Grup_variable <- "Score_pelvic" #### Name of group variable
Analysis_name <- "Fecal_final" ###Group analysis name
Results_directory <- "./Clean_codes/Output_Pelvic"
Normalitzation <- "TMM"


for (analysis in list.files(paste0(getwd(),"/Clean_codes/"), pattern = ".R")) { ####Put the directory whith all functions
  source(paste0(getwd(),"/Clean_codes/",analysis))
}

##########################################
##################Fecal###################
##########################################


############
##### Plot composition function
############


Composition_function(Gut, ####Physeq
                     Grup_variable, ###Group variable
                     "Genus", ###Tax graf
                     Analysis_name, ###File_names
                     Results_directory) ## ## Result of file names


############
##### Especific taxa composition function
############

Organism <- c("Escherichia coli","Enterococcus","Staphylococcus","Klebsiella","Proteus")



for (i in Organism) {
  
  for (Micro in unique(rownames(otu_table(Gut))[grep(pattern = i,rownames(otu_table(Gut)))])) {
    Taxa_function(Physeq = Gut,Grup_variable = Grup_variable,Microorganism = Micro,
                  Analysis_name = Analysis_name,
                  Results_directory = Results_directory
                    )
  }
}


############
##### Alpha diversity function
############

Gut_rarefied = rarefy_even_depth(Gut, rngseed=1, sample.size=0.99*min(sample_sums(Vaginal)), replace=F)

Alpha_Gut <- Alpha_function(Physeq = Gut,Diversity_index = c("Observed","Chao1","Shannon","Simpson","InvSimpson" ,"Fisher"),
               Create_normal_plot = "No",
               Grup_variable = Grup_variable,
               Adjformula = Adjformula,
               Create_diversity_plot = "No")



############
##### Beta diversity function
############
Beta_function(Gut,
              Adjformula,
              "TMM", ###'"relative","TMM","vst","log2".
              "bray") ####### Could be bray, euclidean for low alfa diversity try Weighted Unifrac, could be better having a low rare species



############
##### Differential ANCOMBC function
############


##### "UTI-no-postcoital"

  Feceslfc <- list()
  
  for (i in c("Order" ,  "Family" , "Genus"   ,"Species")) {
    
    
    Feceslfc <- append(Feceslfc,Diferential_function(Gut,Taxonomy = paste(i),
                                                     "                                                             bonferroni",
                                                     Adjformula,
                                                     Grup_variable,
                                                     Results_directory))
  }
  
  Results_list_lfc1 <- list(Feceslfc[[1]],Feceslfc[[3]],Feceslfc[[5]],Feceslfc[[7]])
  
  library(openxlsx)
  
  wb <- createWorkbook()
  
  for (i in c("Order" ,  "Family" , "Genus"   ,"Species")) {
    
    if(i == "Order"){
      j <- 0
    }
    j <- j+1
    addWorksheet(wb, i)
    Model_data <- (Results_list_lfc1[[j]] %>%
                     dplyr::select(taxon, contains("Pelvic")))[,1:6]
    Model_data$BH_Pvalue <- p.adjust(p = Model_data$p_Score_pelvic, method = "BH")
    colnames(Model_data) <- c("Taxon","LFC","SE","W","Pvalue","Bonferroni_Pvalue","BH_Pvalue")
    writeData(wb, sheet = i, x = Model_data)
  }
  saveWorkbook(wb, file = "./Clean_codes/Articles/Fecal_model_score_pelvic.xlsx", overwrite = TRUE)
  
  
  
  
  
  Plot_list_lfc <- list(Feceslfc[[2]],Feceslfc[[4]],Feceslfc[[6]],Feceslfc[[8]])
  
  
  plot_complete <- ggarrange(plotlist = Plot_list_lfc,labels = c("Order" ,  "Family" , "Genus"   ,"Species"),ncol = 2,
                             nrow = 2,common.legend = T)
  
  ggsave(plot = plot_complete,paste0(Results_directory,"/",Analysis_name,"pelvic_score.png"), width = 20, height =20, units = "cm")
  



##########################################
##################Vaginal#################
##########################################

Analysis_name <- "Vaginal_final" ###Group analysis name


############
##### Plot composition function
############


Composition_function(Vaginal, ####Physeq
                     Grup_variable, ###Group variable
                     "Genus", ###Tax graf
                     Analysis_name, ###File_names
                     Results_directory) ## ## Result of file names



############
##### Alpha diversity function
############

  
  Vaginal_rarefied = rarefy_even_depth(Vaginal, rngseed=1, sample.size=0.99*min(sample_sums(Vaginal)), replace=F)

Alpha_vaginal <- Alpha_function(Physeq = Vaginal,Diversity_index = c("Observed","Chao1","Shannon","Simpson","InvSimpson" ,"Fisher"),
                            Create_normal_plot = "No",
                            Grup_variable = Grup_variable,
                            Adjformula = Adjformula,
                            Create_diversity_plot = "No")

############
##### Beta diversity function
############
Beta_function(Vaginal,
              Adjformula,
              "TMM", ###'"relative","TMM","vst","log2".
              "bray") ####### Could be bray, euclidean for low alfa diversity try Weighted Unifrac, could be better having a low rare species



############
##### Differential ANCOMBC function
############


  
  ##### "UTI-postcoital"
  
  Feceslfc <- list()
  
  
  for (i in c("Order" ,  "Family" , "Genus"   ,"Species")) {
    
    
    Feceslfc <- append(Feceslfc,Diferential_function(Vaginal,Taxonomy = paste(i),
                                                     "bonferroni",
                                                     Adjformula,
                                                     Grup_variable,
                                                     Results_directory))
  }
  
  
Results_list_lfc1 <- list(Feceslfc[[1]],Feceslfc[[3]],Feceslfc[[5]],Feceslfc[[7]])

library(openxlsx)

wb <- createWorkbook()

for (i in c("Order" ,  "Family" , "Genus"   ,"Species")) {
  
  if(i == "Order"){
    j <- 0
  }
  j <- j+1
  addWorksheet(wb, i)
  Model_data <- (Results_list_lfc1[[j]] %>%
                   dplyr::select(taxon, contains("Pelvic")))[,1:6]
  Model_data$BH_Pvalue <- p.adjust(p = Model_data$p_Score_pelvic, method = "BH")
  colnames(Model_data) <- c("Taxon","LFC","SE","W","Pvalue","Bonferroni_Pvalue","BH_Pvalue")
  writeData(wb, sheet = i, x = Model_data)
}
saveWorkbook(wb, file = "./Clean_codes/Articles/Vaginal_model_score_pelvic.xlsx", overwrite = TRUE)

  
  plot_complete <- ggarrange(plotlist = Plot_list_lfc,labels = c("Order" ,  "Family" , "Genus"   ,"Species"),ncol = 2,
                             nrow = 2,common.legend = T)
  
  ggsave(plot = plot_complete,paste0(Results_directory,"/",Analysis_name,"Pelvic_score.png"), width = 20, height = 20, units = "cm")
  
  

##########################################
##################Both####################
##########################################

############
##### Graf Both diversity to compare the effect beetween Fecal and vaginal microbiome.
############


Models_alpha_gut <- Alpha_Gut[[2]]
Models_alpha_vaginal <- Alpha_vaginal[[2]]

Models_alpha_gut$Microbiome <- "Gut"
Models_alpha_vaginal$Microbiome <- "Vaginal"

Models_alpha <- rbind.data.frame(Models_alpha_vaginal,Models_alpha_gut)

### prepare data

Models_alpha$beta <- as.numeric(Models_alpha$beta)
Models_alpha$se <- as.numeric(Models_alpha$se)

Models_alpha <- Models_alpha[Models_alpha$name != "(Intercept)",]

library(ggforestplot)

list_plots <- list()
z = 1

for (i in c("Chao1","Shannon","InvSimpson")) {
  
  
  Model_filtred <- Models_alpha[Models_alpha$Index == i,]
  
  # Forestplot
  list_plots[[z]] <- forestplot(
    df = Model_filtred,
    estimate = beta,
    logodds = FALSE,
    colour = Microbiome,
    title = paste0("Beta values from ",i,  " index"),
    xlab = "The alpha estimates against the patients",pvalue = pvalue
  )
  z <- z+1
}


plot_complete <- ggarrange(plotlist = list_plots,labels = NULL,ncol = 3,
                           nrow = 1,common.legend = T)


ggsave(plot = plot_complete,paste0(Results_directory,"/",Analysis_name,"_all_alfa_models.png"), width = 20, height = 20, units = "cm")





############
##### Beta diversity between two groups
############

###Merge fecal and vaginal phyloseq

All_phylo <- merge_phyloseq(Gut, Vaginal)


relab_genera <- microbial::normalize(All_phylo, method = "TMM")


#### Perform the distance matrix
wu = phyloseq::distance(relab_genera, "bray")
wu.m = reshape::melt(as.matrix(wu))
wu.m = wu.m %>%
  filter(as.character(X1) != as.character(X2)) %>%
  mutate_if(is.factor, as.character)

sd = data.frame(sample_data(All_phylo))

sd = sd %>%
  select("colnames.mostres_filt.", "Pacient") %>%
  mutate_if(is.factor,as.character)
colnames(sd) = c("colnames.mostres_filt.", "Var1")
wu.sd = merge(wu.m, sd,by.x="X1" ,by.y = "colnames.mostres_filt.")
colnames(sd) = c("colnames.mostres_filt.", "Var2")
wu.sd = merge(wu.sd, sd,by.x="X2" ,by.y = "colnames.mostres_filt.")

wu.sd <- wu.sd %>%
  dplyr::mutate(Plot_grup = case_when(Var1 == "Yes" & Var2 == "Yes" ~ "Pacient",
                                      Var1 == "No" & Var2 == "No" ~ "Control",
                                      Var1 == "Yes" & Var2 == "No" ~ "Control vs Pacient",
                                      Var1 == "No" & Var2 == "Yes" ~ "Control vs Pacient"
                                      
  )) %>%
  dplyr::filter(Plot_grup != "Control vs Pacient") %>%
  dplyr::mutate(Sample1 = substr(X1, start = 1, stop = nchar(X1)-3))%>%
  dplyr::mutate(Sample2 = substr(X2, start = 1, stop = nchar(X2)-3))%>%
  dplyr::filter(as.character(Sample1) == as.character(Sample2))


wu.data <- merge(wu.sd,data.frame(sample_data(Vaginal)),by.x = "Sample1",by.y= "patient")




###perform the plot



p =  ggplot(wu.data[!wu.data$Post_coital == "Control",], aes(x = Post_coital, y = value,
                       color = Post_coital))+
  theme_bw() +
  geom_boxplot()+
  geom_jitter(alpha = 0.2)+
  ggpubr::stat_compare_means()+
  ylab("bray") +
  xlab("type") + 
  scale_y_continuous(limits = c(0.98,1))+
  scale_color_manual(values = c( "#3C5488FF","#fcba03")) 


p

