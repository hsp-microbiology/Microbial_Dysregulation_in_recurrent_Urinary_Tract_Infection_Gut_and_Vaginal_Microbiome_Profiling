###################################################################
###################################################################
###################################################################
###     Functions Main file exam -> Created by Adrià         ######
###  Linkedin: https://www.linkedin.com/in/adria-cruells/    ######
###  ORCID: https://orcid.org/0000-0002-1179-7997            ######
###  Github:    https://github.com/ADRS94                    ######
###################################################################
###################################################################

r

setwd("") ###### put the directory were you work
Vaginal <- readRDS(file = "")  #Enter directory of a created phytloseq( *:RData)
Gut <- readRDS(file = "")  #Enter directory of a created phytloseq( *:RData)    subset_taxa(!Genus == "Homo")%>% for vaginal and fecal microbiome

##########
####Filter
##########
library(dplyr)
library(phyloseq)

Vaginal <- Vaginal %>%
  subset_taxa(Genus != "Homo")%>%
  subset_taxa(Genus != "Ralstonia")%>%
  subset_samples(Pacient=="Yes")

Gut <- Gut %>%
  subset_taxa(!grepl("unk",Genus)) %>%
  subset_samples(Pacient=="Yes")



################
##Analysis parameters
################

Adjformula <- ("Score_pelvic+Age+BMI+Menopause") #####Model formula (Ex model 3)
Grup_variable <- "Score_pelvic" #### Name of group variable
Analysis_name <- "Score_pelvic_fecal" ###Group analysis name
Results_directory <- "./Clean_codes/Output_try_dis"
Normalitzation <- "TMM"


for (analysis in list.files(paste0(getwd(),"/Clean_codes/"), pattern = ".R")) { ####Put the directory whith all functions
  source(paste0(getwd(),"/Clean_codes/",analysis))
}

##########################################
##################Fecal###################
##########################################



############
##### Alpha diversity function
############






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


Feceslfc <- list()
j <- 0

for (i in c("Order" ,  "Family" , "Genus"   ,"Species")) {
  
  
  Feceslfc <- append(Feceslfc,Diferential_function(Gut,Taxonomy = paste(i),
                                                   "fdr",
                                                   Adjformula,
                                                   Grup_variable,
                                                   Results_directory))
}


Plot_list_lfc <- list(Feceslfc[[2]],Feceslfc[[4]],Feceslfc[[6]],Feceslfc[[8]])


plot_complete <- ggarrange(plotlist = Plot_list_lfc,labels = c("Order" ,  "Family" , "Genus"   ,"Species"),ncol = 2,
                           nrow = 2,common.legend = T)


##########################################
##################Vaginal#################
##########################################

Analysis_name <- "Score_pelvic_vaginal" ###Group analysis name


############
##### Alpha diversity function
############


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


vagilfc <- list()
j <- 0

for (i in c("Order" ,  "Family" , "Genus"   ,"Species")) {
  
  
  vagilfc <- append(Feceslfc,Diferential_function(Vaginal,Taxonomy = paste(i),
                                                   "fdr",
                                                   Adjformula,
                                                   Grup_variable,
                                                   Results_directory))
}


Plot_list_lfc <- list(vagilfc[[2]],vagilfc[[4]],vagilfc[[6]],vagilfc[[8]])


plot_complete <- ggarrange(plotlist = Plot_list_lfc,labels = c("Order" ,  "Family" , "Genus"   ,"Species"),ncol = 2,
                           nrow = 2,common.legend = T)


##########################################
##################Both####################
##########################################

############
##### Effect of Gut diversity in vaginal diversity 
############


Alpha_Gut$sample <- substr(Alpha_Gut$Type,1,nchar(Alpha_Gut$Type)-3)
Alpha_vaginal$sample <- substr(Alpha_vaginal$Type,1,nchar(Alpha_vaginal$Type)-3)

colnames(Alpha_Gut) <- paste0(colnames(Alpha_Gut),"_fecal")
colnames(Alpha_vaginal) <- paste0(colnames(Alpha_vaginal),"_vaginal")




Alpha_data <- base::merge(Alpha_Gut,Alpha_vaginal, by.x = "sample_fecal", by.y = "sample_vaginal")




model <- glm(Shannon_vaginal ~ Shannon_fecal
             , data = Alpha_data)



model <- glm(Shannon_fecal ~Shannon_vaginal 
             , data = Alpha_data)


summary(model)


plot(model)
