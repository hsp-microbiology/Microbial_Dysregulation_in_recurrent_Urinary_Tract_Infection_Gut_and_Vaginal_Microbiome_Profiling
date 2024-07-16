###################################################################
###################################################################
###################################################################
###     Functions Main file exam -> Created by Adrià         ######
###  Linkedin: https://www.linkedin.com/in/adria-cruells/    ######
###  ORCID: https://orcid.org/0000-0002-1179-7997            ######
###  Github:   https://github.com/ADRS94                     ######
###################################################################
###################################################################


######## Enter the phyloseq object, set working directory and chargue libraries


setwd("") ###### put the directory were you work
Vaginal <- readRDS(file = "")  #Enter directory of a created phytloseq( *:RData)
Gut <- readRDS(file = "")  #Enter directory of a created phytloseq( *:RData)  

##########
####Filter the non-classified species
##########
library(dplyr)

Vaginal <- Vaginal %>%
  subset_taxa(!grepl("unk",Genus))%>%
  subset_taxa(!grepl("unk",Species))

Gut <- Gut %>%
  subset_taxa(!grepl("unk",Genus))%>%
  subset_taxa(!grepl("unk",Species))



##########################
####Patients data analysis
##########################
library(table1)

#### obtain the data from phyloseq and select the interest variables ( in this case (Patient//Age(years)//BMI//Menopause//UTI/year))

Patients_data <- data.frame(sample_data(Vaginal)) %>%
  dplyr::select(Pacient,BMI,Age,Menopause,Score_pelvic,Postcoital.UTI) 

# Patients_data$Score_pelvic <- as.numeric(Patients_data$Score_pelvic)
# Patients_data$Score_pelvic[Patients_data$Pacient == "No"] <- 0

Patients_data$Postcoital.UTI <- as.factor(Patients_data$Postcoital.UTI)
Patients_data$Postcoital.UTI[Patients_data$Pacient == "No"] <- "No"


# Factor the basic variables that
# we're interested in

####Rename all the variables to get the correct table 

Patients_data$Menopause <- 
  factor(Patients_data$Menopause, levels=c("Yes",
                                           "No"),
         labels=c("Yes", 
                  "No"))

label(Patients_data$BMI)       <- "Body mass index (BMI)"
label(Patients_data$Menopause)       <- "Menopause"
label(Patients_data$Age)     <- "Age (Years)"
label(Patients_data$Score_pelvic) <- "Score pelvic"
label(Patients_data$Postcoital.UTI) <- "UTI postcoital"



table1(~  BMI + Menopause + Age + Postcoital.UTI + Score_pelvic   | Pacient, data=Patients_data)

################
##Analysis parameters
################

Adjformula <- ("Pacient+Age+BMI+Menopause") #####Model formula 
Grup_variable <- "Pacient" #### Name of group variable
Analysis_name <- "Fecal" ###Group analysis name
Results_directory <- ""
Normalitzation <- "TMM" ## For beta diversity


for (analysis in list.files(paste0(getwd(),"/Clean_codes/"), pattern = ".R")) { ####Put the directory whith all functions
  source(paste0(getwd(),"/Clean_codes/",analysis))
}
library(microViz)

## Modify metadata names

Gut <- Gut %>%
  ps_mutate(
    Pacient = case_when(Pacient == "No" ~ "wctrl",
                            Pacient == "Yes" ~ "wrUTI")
  )

Vaginal <- Vaginal %>%
  ps_mutate(
    Pacient = case_when(Pacient == "No" ~ "wctrl",
                        Pacient == "Yes" ~ "wrUTI")
  )




##########################################
##################Fecal###################
##########################################


############
##### Plot composition function
############


Fig1_A <- Composition_function(Gut, ####Physeq
                     Grup_variable, ###Group variable
                     "Genus", ###Tax graf
                     Analysis_name, ###File_names
                     Results_directory) ## ## Result of file names



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

Fig1_D <- Beta_function(Gut,
              Adjformula,
              "TMM", ###'"relative","TMM","vst","log2".
              "bray",
              Grup_variable = Grup_variable) ####### Could be bray, euclidean for low alfa diversity try Weighted Unifrac, could be better having a low rare species



############
##### Differential ANCOMBC function
############
Analysis_name <- "Fecal" ###Group analysis name

Feceslfc <- list()
j <- 0

for (i in c("Order" ,  "Family" , "Genus"   ,"Species")) {
  
  
  Feceslfc <- append(Feceslfc,Diferential_function(Gut,Taxonomy = paste(i),
                                                   "bonferroni",
                                                   Adjformula,
                                                   Grup_variable,
                                                   Results_directory))
  
  
  
  
}


Plot_list_lfc <- list(Feceslfc[[2]],Feceslfc[[4]],Feceslfc[[6]],Feceslfc[[8]])
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
    dplyr::select(taxon, contains(Grup_variable)))[,1:6]
  Model_data$BH_Pvalue <- p.adjust(p = Model_data$p_PacientYes, method = "BH")
  colnames(Model_data) <- c("Taxon","LFC","SE","W","Pvalue","Bonferroni_Pvalue","BH_Pvalue")
  writeData(wb, sheet = i, x = Model_data)
}


saveWorkbook(wb, file = "./Clean_codes/Articles/Fecal_model.xlsx", overwrite = TRUE)


plot_complete <- ggarrange(plotlist = Plot_list_lfc,labels = c("Order" ,  "Family" , "Genus"   ,"Species"),ncol = 2,
                           nrow = 2,common.legend = T)


ggsave(plot = plot_complete,paste0(Results_directory,"/",Analysis_name,"_all_BH_models_Bonferroni.png"), width = 30, height = 30, units = "cm")



#############
##### Discriminant analysis ( Biomarker search) coda 4 microbiome Not function outside, error runing the model
#############

library(coda4microbiome)


set.seed(123) # to reproduce the results

# Get matric from phyloseq

Gut_Genus <-tax_agg(Gut,rank = "Genus") ##from microviz package


Gut_matrix <- data.frame(t(data.frame(otu_table(Gut_Genus)))) ### Get genus possibly, check for the better AUC
Gut_sample <- factor(sample_data(Gut_Genus)$Pacient)
Gut_covar <- data.frame(sample_data(Gut_Genus)$Age,sample_data(Gut_Genus)$BMI,
                        factor(sample_data(Gut_Genus)$Menopause))
Gut_sample <- relevel(Gut_sample, ref = "No")

coda_glmnet_Crohn<-coda_glmnet(x=Gut_matrix,y=Gut_sample,covar = Gut_covar,nfolds = 5) ## Care sort the mattrix as sample_data 

saveRDS(coda_glmnet_Crohn,"./Coda_Gut_species.rds")



coda_glmnet_Crohn$taxa.num

coda_glmnet_Crohn$taxa.name

coda_glmnet_Crohn$`log-contrast coefficients`

sum(coda_glmnet_Crohn$`log-contrast coefficients`) ###Should be 0 


coef<-coda_glmnet_Crohn$`log-contrast coefficients`
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
coda_glmnet_Crohn$taxa.name[positives[op]]


negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
coda_glmnet_Crohn$taxa.name[negatives[on]]


coda_glmnet_Crohn$`signature plot`


##########################################
##################Vaginal#################
##########################################

Analysis_name <- "Vaginal" ###Group analysis name


############
##### Plot composition function
############


Fig1_B <- Composition_function(Vaginal, ####Physeq
                     Grup_variable, ###Group variable
                     "Genus", ###Tax graf
                     Analysis_name, ###File_names
                     Results_directory) ## ## Result of file names



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
Fig1_E <- Beta_function(Vaginal,
              Adjformula,
              "TMM", ###'"relative","TMM","vst","log2".
              "bray") ####### Could be bray, euclidean for low alfa diversity try Weighted Unifrac, could be better having a low rare species



############
##### Differential ANCOMBC function
############


vagilfc <- list()
j <- 0

for (i in c("Order" ,  "Family" , "Genus"   ,"Species")) {
  
  
  vagilfc <- append(vagilfc,Diferential_function(Vaginal,Taxonomy = paste(i),
                                                   "bonferroni",
                                                   Adjformula,
                                                   Grup_variable,
                                                   Results_directory))
}


Results_list_lfc1 <- list(vagilfc[[1]],vagilfc[[3]],vagilfc[[5]],vagilfc[[7]])

library(openxlsx)

wb <- createWorkbook()

for (i in c("Order" ,  "Family" , "Genus"   ,"Species")) {
  
  if(i == "Order"){
    j <- 0
  }
  j <- j+1
  addWorksheet(wb, i)
  Model_data <- (Results_list_lfc1[[j]] %>%
                   dplyr::select(taxon, contains(Grup_variable)))[,1:6]
  Model_data$BH_Pvalue <- p.adjust(p = Model_data$p_PacientYes, method = "BH")
  colnames(Model_data) <- c("Taxon","LFC","SE","W","Pvalue","Bonferroni_Pvalue","BH_Pvalue")
  writeData(wb, sheet = i, x = Model_data)
}
saveWorkbook(wb, file = "./Clean_codes/Articles/Vaginal_model.xlsx", overwrite = TRUE)



Plot_list_lfc <- list(vagilfc[[2]],vagilfc[[4]],vagilfc[[6]],vagilfc[[8]])


plot_complete <- ggarrange(plotlist = Plot_list_lfc,labels = c("Order" ,  "Family" , "Genus"   ,"Species"),ncol = 2,
                           nrow = 2,common.legend = T)


ggsave(plot = plot_complete,paste0(Results_directory,"/",Analysis_name,"_all_bonferroni_models.png"), width = 30, height = 30, units = "cm")



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

for (i in c("Chao1","Shannon","Simpson")) {
  
  
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


ggsave(plot = plot_complete,paste0(Results_directory,"/",Analysis_name,"_all_alfa_models.png"), width = 40, height = 20, units = "cm")


############
##### Graf Both diversity to compare the effect beetween Fecal and vaginal microbiome.
############


Models_alpha_gut <- Alpha_Gut_rarefied[[2]]
Models_alpha_vaginal <- Alpha_vaginal_rarefied[[2]]

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

for (i in c("Chao1","Shannon","Simpson")) {
  
  
  Model_filtred <- Models_alpha[Models_alpha$Index == i,]
  
  # Forestplot
  list_plots[[z]] <- forestplot(
    df = Model_filtred,
    estimate = beta,
    logodds = FALSE,
    colour = Microbiome,
    title = paste0(i),
    xlab = "",pvalue = pvalue
  )
  z <- z+1
}


plot_complete <- ggarrange(plotlist = list_plots,labels = NULL,ncol = 3,
                           nrow = 1,common.legend = T,legend = "none")


ggsave(plot = plot_complete,paste0(Results_directory,"/",Analysis_name,"_rarefy_alfa_models.png"), width = 40, height = 20, units = "cm")



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
  dplyr::mutate(Plot_grup = case_when(Var1 == "wrUTI" & Var2 == "wrUTI" ~ "wrUTI",
                                      Var1 == "wctrl" & Var2 == "wctrl" ~ "wctrl",
                                      Var1 == "wrUTI" & Var2 == "wctrl" ~ "wctrl vs wrUTI",
                                      Var1 == "wctrl" & Var2 == "wrUTI" ~ "wctrl vs wrUTI"
    
  )) %>%
  dplyr::filter(Plot_grup != "wctrl vs wrUTI") %>%
  dplyr::mutate(Sample1 = substr(X1, start = 1, stop = nchar(X1)-3))%>%
  dplyr::mutate(Sample2 = substr(X2, start = 1, stop = nchar(X2)-3))%>%
  dplyr::filter(as.character(Sample1) == as.character(Sample2))





###perform the plot



Fig1_F =  ggplot(wu.sd, aes(x = Plot_grup, y = value,
                       color = Plot_grup))+
  theme_bw() +
  geom_boxplot()+
  geom_jitter(alpha = 0.2)+
  #ggpubr::stat_compare_means()+
  ylab("bray") +
  xlab("type") + 
  scale_y_continuous(limits = c(0.98,1))+
  scale_color_manual(values = c("#3C5488FF","#00A087FF")) 

p


ggsave(plot = p,paste0(Results_directory,"/",Analysis_name,"Beta_vaginal_fecal.png"), width = 40, height = 20, units = "cm")


######Selbal for microbiome in fecal 

library(selbal)

# Define x, y and z
All_dat <- merge(data.frame(t(otu_table(Gut %>%
                                          tax_agg("Genus")))),data.frame(sample_data(Gut)),by = "row.names" ) %>%
  column_to_rownames(var = "Row.names")



noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- colSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[,bigones]
  print(percent)
  return(Matrix_1)
}




x <- as.data.frame(noise.removal(All_dat[,1:1052])) ###841 for genus
y <- as.factor(All_dat[,1126]) ## Group study
z <- data.frame(Age = All_dat[,1056],
                BMI = All_dat[,1060],
                Menopause =All_dat[,1081]
                #Sexual_frequency =All_dat[,922],
                #Vaginal_Birth =All_dat[,921]
                ) ## Covariables
CV.BAL.dic <- selbal.cv(x = x, y = y, n.fold = 5, n.iter = 10,covar = z,
                       logit.acc = "AUC") ######The columns with more than 80% of zeros gimes error.



CV.BAL.dic$accuracy.nvar ### The optimal is 2,



Fig2_A<- CV.BAL.dic$global.plot


####Vaginal



# Define x, y and z
All_dat <- merge(data.frame(t(otu_table(Vaginal %>%
                                          tax_agg("Genus")))),data.frame(sample_data(Vaginal)),by = "row.names" ) %>%
  column_to_rownames(var = "Row.names")



noise.removal <- function(dataframe, percent=0.1, top=NULL){
  dataframe->Matrix
  bigones <- colSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[,bigones]
  print(percent)
  return(Matrix_1)
}




x <- as.data.frame(noise.removal(All_dat[,1:486])) ###841 for genus
y <- as.factor(All_dat[,560]) ## Group study
z <- data.frame(Age = All_dat[,490],
                BMI = All_dat[,494],
                Menopause =All_dat[,515]
                #Sexual_frequency =All_dat[,922],
                #Vaginal_Birth =All_dat[,921]
) ## Covariables
CV.BAL.val <- selbal.cv(x = x, y = y, n.fold = 5, n.iter = 10,covar = z,
                        logit.acc = "AUC") ######The columns with more than 80% of zeros gimes error.



CV.BAL.val$accuracy.nvar ### The optimal is 2,


Fig2_B<- CV.BAL.val$global.plot




### Figure 1 

# Show the figure 1

Fig1 <- ggarrange(
  ggarrange(Fig1_A, Fig1_B, ncol = 2,
            labels = c("A", "B")),
  plot_complete,
  ggarrange(Fig1_D, Fig1_E, ncol = 2,
            labels = c("D", "E")),
  Fig1_F,
  nrow = 4,
  labels = c("", "C", "", "F")
)

Fig1 <- annotate_figure(Fig1,
                                 top = text_grob("Figure 1", color = "black", face = "bold", size = 14))


# Guardar la figura como una imagen de alta calidad
ggsave("./../Article_ultima_revisió/Figures/Figure1.png", plot = Fig1, width = 12, height = 14, dpi = 300) ###DPI mayor, mayor calidad


#### Show the figure 2

Fig2 <- ggarrange(Fig2_A, Fig2_B, ncol = 2,
          labels = c("A", "B"))

Fig2 <- annotate_figure(Fig2,
                        top = text_grob("Figure 2", color = "black", face = "bold", size = 14))


# Guardar la figura como una imagen de alta calidad
ggsave("./../Article_ultima_revisió/Figures/Figure2.png", plot = Fig2, width = 18, height = 10, dpi = 300) ###DPI mayor, mayor calidad



#### Show the figure 3

### Read the tables:

library(readxl)




Table_gut_1 <- read_excel(path = "./../Article_ultima_revisió/Supplementary tables/Supplementary_tables_1.xlsx",sheet = "DA_species_gut",skip = 2)



Table_gut_1 <-  Table_gut_1 %>%
  dplyr::mutate(diffexpressed   = case_when(`P-value` >= 0.05 ~ "No",
                                            `P-value` <= 0.05 ~ "Yes"
  ),
  delabel = case_when(`Adj. P-value (BH)` >= 0.05 ~ NA,
                      `Adj. P-value (BH)` <= 0.05 ~ Species
  )
  )


Fig3_A <- ggplot(data=Table_gut_1, aes(x=Effect, y=-log10(`P-value`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept= 0, col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")





Table_gut_1 <- read_excel(path = "./../Article_ultima_revisió/Supplementary tables/Supplementary_tables_1.xlsx",sheet = "DA_genus_gut",skip = 2)



Table_gut_1 <-  Table_gut_1 %>%
  dplyr::mutate(diffexpressed   = case_when(`P-value` >= 0.05 ~ "No",
                                            `P-value` <= 0.05 ~ "Yes"
  ),
  delabel = case_when(`Adj. P-value (BH)` >= 0.05 ~ NA,
                      `Adj. P-value (BH)` <= 0.05 ~ Genera
  )
  )


Fig3_B <- ggplot(data=Table_gut_1, aes(x=Effect, y=-log10(`P-value`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept= 0, col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")




Table_gut_1 <- read_excel(path = "./../Article_ultima_revisió/Supplementary tables/Supplementary_tables_1.xlsx",sheet = "DA_species_vaginal",skip = 2)



Table_gut_1 <-  Table_gut_1 %>%
  dplyr::mutate(diffexpressed   = case_when(`P-value` >= 0.05 ~ "No",
                                            `P-value` <= 0.05 ~ "Yes"
  ),
  delabel = case_when(`Adj. P-value (BH)` >= 0.05 ~ NA,
                      `Adj. P-value (BH)` <= 0.05 ~ Species
  )
  )


Fig3_C <- ggplot(data=Table_gut_1, aes(x=Effect, y=-log10(`P-value`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept= 0, col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")




Table_gut_1 <- read_excel(path = "./../Article_ultima_revisió/Supplementary tables/Supplementary_tables_1.xlsx",sheet = "DA_genus_vaginal",skip = 2)



Table_gut_1 <-  Table_gut_1 %>%
  dplyr::mutate(diffexpressed   = case_when(`P-value` >= 0.05 ~ "No",
                                            `P-value` <= 0.05 ~ "Yes"
  ),
  delabel = case_when(`Adj. P-value (BH)` >= 0.05 ~ NA,
                      `Adj. P-value (BH)` <= 0.05 ~ Genera
  )
  )


Fig3_D <- ggplot(data=Table_gut_1, aes(x=Effect, y=-log10(`P-value`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept= 0, col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")





Fig3 <- ggarrange(
  ggarrange(Fig3_A, Fig3_B, ncol = 2,
            labels = c("A", "B"),common.legend = T),
  ggarrange(Fig3_C, Fig3_D, ncol = 2,
            labels = c("C", "D"),common.legend = T),
  nrow = 2,common.legend = T
)

Fig3 <- annotate_figure(Fig3,
                        top = text_grob("Figure 3", color = "black", face = "bold", size = 14))


# Guardar la figura como una imagen de alta calidad
ggsave("./../Article_ultima_revisió/Figures/Figure3.png", plot = Fig3, width = 18, height = 14, dpi = 300) ###DPI mayor, mayor calidad
