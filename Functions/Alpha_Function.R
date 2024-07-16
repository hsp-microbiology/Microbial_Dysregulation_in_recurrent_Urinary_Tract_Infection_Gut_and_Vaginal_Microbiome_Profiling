###################################################################
###################################################################
###################################################################
###            Functions ANCOMBC -> Created by Adrià         ######
###  Linkedin: https://www.linkedin.com/in/adria-cruells/    ######
###  ORCID: https://orcid.org/0000-0002-1179-7997            ######
###  Github:    Working on it                                ######
###################################################################
###################################################################

###################################################
###############Chargue libraries###################
###################################################


pacman::p_load(ggpubr,####Analysis package
               caret,#### Analysis package
               phyloseq,#### Filetype package
               readxl, #### Entry table package
               openxlsx, #### Entry table package
               dplyr, #### Table manage package
               tidyverse, #### Plot package
               ggrepel, #### Plot package
               sjPlot #####Plot package
               ) 




###################################################
###############Alpha   FUNCTIONS###################
###################################################


Alpha_function <- function(Physeq,Diversity_index,Create_normal_plot,Grup_variable,Adjformula,Create_diversity_plot){
  
  
  ############################
  ### Alpha computing run  ###
  ############################

  
  data_alpha <- estimate_richness(Physeq, split = TRUE, measures = NULL)
  data_alpha$Type <- rownames(data_alpha)
  Alpha_computed <- merge(data_alpha,as.data.frame(sample_data(Physeq)), by = "row.names")

  ##################
  ##Create plot  ###
  ##################

  if(Create_normal_plot == "Yes"){
    for (i in Diversity_index) {
      p1 <- ggviolin(Alpha_computed, x = Grup_variable, y = i,
                     add = "boxplot", fill = Grup_variable, palette = c("#00A087FF", "#3C5488FF","#fcba03")) 
      p1 <- p1 + stat_compare_means() 
      ggsave(plot = p1,paste0(Results_directory,"/",Analysis_name,i,"raw.png"), width = 20, height = 20, units = "cm")
      

    }
    
  }

  
  ##################################################
  ##Run the selected models and plot the models  ###
  ##################################################

  ###############
  ##Observed  ###
  ###############
  
  if("Observed" %in% Diversity_index){
  
    ModelEquation <- as.formula(paste0("Observed ~",Adjformula))
  
  
  model_Observed <- glm(ModelEquation
               , data = Alpha_computed)
  
  sink(file = paste0(Results_directory,"/",Analysis_name,"Observed.txt"))
  print(summary(model_Observed))
  sink()
  Plot_Observed <-  plot_model(model_Observed)
  png(filename = paste0(Results_directory,"/",Analysis_name,"Observed.png"),
      width = 800, height = 800, units = "px")
  print(plot_model(model_Observed))
   dev.off()
  
  }
  
  ############
  ##Chao1  ###
  ############
  
  if("Chao1" %in% Diversity_index){
   
     ModelEquation <- as.formula(paste0("Chao1 ~",Adjformula))
    
    
    model_Chao1 <- glm(ModelEquation
                 , data = Alpha_computed)
    
    sink(file = paste0(Results_directory,"/",Analysis_name,"Chao1.txt"))
    print(summary(model_Chao1))
    sink()
    Plot_Chao1 <-  plot_model(model_Chao1)
    png(filename = paste0(Results_directory,"/",Analysis_name,"Chao1.png"),
        width = 800, height = 800, units = "px")
    print(plot_model(model_Chao1))
    dev.off()
    
  }
  
  ##########
  ##ACE  ###
  ##########
  
  
  if("ACE" %in% Diversity_index){
    
    ModelEquation <- as.formula(paste0("ACE ~",Adjformula))

    model <- glm(ModelEquation
                 , data = Alpha_computed)
    
    sink(file = paste0(Results_directory,"/",Analysis_name,"ACE.txt"))
    print(summary(model))
    sink()
    Plot_ACE <-  plot_model(model)
    png(filename = paste0(Results_directory,"/",Analysis_name,"ACE.png"),
        width = 800, height = 800, units = "px")
    print(plot_model(model))
    dev.off()
    
  }
  ##############
  ##Shannon  ###
  ##############
  
  if("Shannon" %in% Diversity_index){
    
    ModelEquation <- as.formula(paste0("Shannon ~",Adjformula))
    
    
    model_Shannon <- glm(ModelEquation
                 , data = Alpha_computed)
    
    sink(file = paste0(Results_directory,"/",Analysis_name,"Shannon.txt"))
    print(summary(model_Shannon))
    sink()
    Plot_Shannon <-  plot_model(model_Shannon)
    png(filename = paste0(Results_directory,"/",Analysis_name,"Shannon.png"),
        width = 800, height = 800, units = "px")
    print(plot_model(model_Shannon))
    dev.off()
    
  }
  
  ##############
  ##Simpson  ###
  ##############
  
  if("Simpson" %in% Diversity_index){
    
    
    ModelEquation <- as.formula(paste0("Simpson ~",Adjformula))

    model_Simpson <- glm(ModelEquation
                 , data = Alpha_computed)
    
    sink(file = paste0(Results_directory,"/",Analysis_name,"Simpson.txt"))
    print(summary(model_Simpson))
    sink()
    Plot_Simpson <-  plot_model(model_Simpson)
    png(filename = paste0(Results_directory,"/",Analysis_name,"Simpson.png"),
        width = 800, height = 800, units = "px")
    print(plot_model(model_Simpson))
    dev.off()
    
  }
  
  #######################
  ##InvSimpson  #########
  #######################
  
  if("InvSimpson" %in% Diversity_index){
    
    
    ModelEquation <- as.formula(paste0("InvSimpson ~",Adjformula))
    
    model_InvSimpson <- glm(ModelEquation
                 , data = Alpha_computed)
    
    sink(file = paste0(Results_directory,"/",Analysis_name,"InvSimpson.txt"))
    print(summary(model_InvSimpson))
    sink()
    Plot_InvSimpson <-  plot_model(model_InvSimpson)
    png(filename = paste0(Results_directory,"/",Analysis_name,"InvSimpson.png"),
        width = 800, height = 800, units = "px")
    print(plot_model(model_InvSimpson))
    dev.off()
    
  }
  
  
  ##########################
  #######Fisher  ###########
  ##########################
  
  if("Fisher" %in% Diversity_index){
    
    
    ModelEquation <- as.formula(paste0("Fisher ~",Adjformula))
    
    model_Fisher <- glm(ModelEquation
                 , data = Alpha_computed)
    
    sink(file = paste0(Results_directory,"/",Analysis_name,"Fisher.txt"))
    print(summary(model_Fisher))
    sink()
    
   Plot_Fisher <-  plot_model(model_Fisher)
    
    png(filename = paste0(Results_directory,"/",Analysis_name,"Fisher.png"),
        width = 800, height = 800, units = "px")
    print( plot_model(model_Fisher))
    dev.off()
    
    
    
    Fisher_estimates <- cbind.data.frame("Fisher",summary(model_Fisher)$coef,rownames(summary(model_Fisher)$coef))
    Shannon_estimates <- cbind.data.frame("Shannon",summary(model_Shannon)$coef,rownames(summary(model_Shannon)$coef))
    InvSimpson_estimates <- cbind.data.frame("InvSimpson",summary(model_InvSimpson)$coef,rownames(summary(model_InvSimpson)$coef))
    Observed_estimates <- cbind.data.frame("Observed",summary(model_Observed)$coef,rownames(summary(model_Observed)$coef))
    Simpson_estimates <- cbind.data.frame("Simpson",summary(model_Simpson)$coef,rownames(summary(model_Simpson)$coef))
    Chao1_estimates <- cbind.data.frame("Chao1",summary(model_Chao1)$coef,rownames(summary(model_Chao1)$coef))
    
    colnames(Fisher_estimates) <- c("Index","beta","se","t_value","pvalue","name")
    colnames(Shannon_estimates) <- c("Index","beta","se","t_value","pvalue","name")
    colnames(InvSimpson_estimates) <- c("Index","beta","se","t_value","pvalue","name")
    colnames(Observed_estimates) <- c("Index","beta","se","t_value","pvalue","name")
    colnames(Simpson_estimates) <- c("Index","beta","se","t_value","pvalue","name")
    colnames(Chao1_estimates) <- c("Index","beta","se","t_value","pvalue","name")  
    Data_forest <- as.tibble(rbind.data.frame(Fisher_estimates,Shannon_estimates,InvSimpson_estimates,Observed_estimates,Simpson_estimates,Chao1_estimates))
    
  }
  

  if(Create_diversity_plot == "Yes"){
    
    All_plots <- ls()[grep("Plot_",ls())]
    All_list_plots <- list()
    ### El error esta en el append
    All_list_plots[[1]] <- Plot_Chao1
    All_list_plots[[2]] <- Plot_Fisher
    All_list_plots[[3]] <- Plot_InvSimpson
    All_list_plots[[4]] <- Plot_Observed
    All_list_plots[[5]] <- Plot_Shannon
    All_list_plots[[6]] <- Plot_Simpson
    plot_complete <- ggarrange(plotlist = All_list_plots,labels = NULL,ncol = 3,
              nrow = 2)
    

    ggsave(plot = plot_complete,paste0(Results_directory,"/",Analysis_name,"_all_models.png"), width = 20, height = 20, units = "cm")
    
    
   
    
    

    

    
    #### make the forest plot
    library(tidyverse)
    # devtools::install_github("NightingaleHealth/ggforestplot")
    library(ggforestplot)
    
    
    Data_forest$beta <- as.numeric(Data_forest$beta)
    Data_forest$se <- as.numeric(Data_forest$se)
    
    Data_forest <- Data_forest[Data_forest$name != "(Intercept)",]
    
#https://nightingalehealth.github.io/ggforestplot/index.html

    # Forestplot
forestplot(
  df = Data_forest,
  estimate = beta,
  logodds = FALSE,
  colour = Index,
  title = "Beta values from different index",
  xlab = "The alpha estimates against the patients",pvalue = pvalue
)

 }
  
  
  
  

  return(list(data_alpha,Data_forest))
  
}



###################################################
####RUN Alpha pipeline(For non source users)#######
###################################################


# Alpha_function(Physeq,Diversity_index,Create_normal_plot,Grup_variable,Adjformula)
