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
Vaginal <- readRDS(file = "C:/Users/Adrià/Desktop/Trabajo/Microbioma R/Clean_codes/Vaginal2.rds")  #Enter directory of a created phytloseq( *:RData)
Gut <- readRDS(file = "C:/Users/Adrià/Desktop/Trabajo/Microbioma R/Clean_codes/Gut.rds")  #Enter directory of a created phytloseq( *:RData)    subset_taxa(!Genus == "Homo")%>% for vaginal and fecal microbiome

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


# > table(sample_data(Vaginal)$Score_pelvic)
# 
# 0  1  2  3  4  5  6  7  8 10 11 12 13 14 15 16 17 18 19 21 22 24 25 27 28 29 31 32 38 39 44 
# 1  1  1  1  4  4  2  2  2  2  1  1  4  3  1  1  2  1  1  1  1  1  1  1  1  1  1  1  1  1  1 
####Vamos a generar dos grupos, aquellas muestras con una puntuacion menor a 10 i con puntuacion major a 10 (20 menos, 27 mas)


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
##### Especific taxa composition function
############

# Organism <- c("Escherichia coli","Staphylococcus","Klebsiella","Proteus")
# 
# 
# All_models <- data.frame()
# 
# 
# 
# for (i in Organism) {
#   
#   for (Micro in unique(rownames(otu_table(Gut))[grep(pattern = i,rownames(otu_table(Gut)))])) {
#    Model_linear <-  Taxa_function_lineal(Physeq = Gut,Adjformula = Adjformula,Microorganism = Micro,
#                   Analysis_name = Analysis_name,
#                   Results_directory = Results_directory
#                     )
#    
#    Linear_effect <- as.data.frame(coef(summary(Model_linear)))
#    
#    All_models <- rbind(All_models,cbind(paste(strsplit(Micro, ";")[[1]][7]),Linear_effect[rownames(Linear_effect) == Grup_variable,]))
#    
#    
#   }
# }
# 
# ###Make the forest plot
# colnames(All_models) <- c("Organism","Estimate","SE","Z_value","P_value")
# write.csv(All_models,file = paste0(Results_directory,"/All_Models",Analysis_name,".csv"))


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

Analysis_name <- "First_try_Vaginal" ###Group analysis name


############
##### Plot composition function
############


# Composition_function(Vaginal, ####Physeq
#                      Grup_variable, ###Group variable
#                      "Genus", ###Tax graf
#                      Analysis_name, ###File_names
#                      Results_directory) ## ## Result of file names


############
##### Especific taxa composition function
############

# Organism <- c("Escherichia coli","Enterococcus","Staphylococcus","Klebsiella","Proteus")
# 
# 
# 
# for (i in Organism) {
#   
#   for (Micro in unique(rownames(otu_table(Vaginal))[grep(pattern = i,rownames(otu_table(Vaginal)))])) {
#     Taxa_function(Physeq = Vaginal,Grup_variable = Grup_variable,Microorganism = Micro,
#                   Analysis_name = Analysis_name,
#                   Results_directory = Results_directory
#     )
#   }
# }


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

############
##### Differential Secom function
############

##### Select the diferential taxa and observe the effect in the other biome.
sigtaxa_vagi <- paste0("Vaginal - ",vagilfc$taxon[vagilfc$diff_PacientYes == "TRUE"])
sigtaxa_feces <- paste0("Fecal - ",Feceslfc$taxon[Feceslfc$diff_PacientYes == "TRUE"])

sigtaxa <- c(sigtaxa_vagi,sigtaxa_feces)
#####Obtain the correlation matrix (the retunr of Both_microbiomes function)

cor_matrix <- as.data.frame(Differential_ecosistem(Gut,Vaginal,"Genus"))

corr_sig <- cor_matrix[grep(paste(sigtaxa,collapse="|"),rownames(cor_matrix)),]


################################ Obtain the correlation only of the Ancombc taxa



