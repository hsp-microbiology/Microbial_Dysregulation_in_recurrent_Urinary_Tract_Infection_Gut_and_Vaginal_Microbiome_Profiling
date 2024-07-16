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
###############ANCOMBC FUNCTIONS###################
###################################################


Diferential_function <- function(Physeq,Taxonomy,Padj_Method,Adjformula,Grup_variable,Results_directory){
  
  
  ##################
  ### Model run  ###
  ##################
  
  
  output = ancombc2(data = Physeq, assay_name = "counts", tax_level = Taxonomy,
                    fix_formula = Adjformula, rand_formula = NULL,
                    p_adj_method = Padj_Method, pseudo = 0, pseudo_sens = FALSE,
                    prv_cut = 0.10, lib_cut = 0.01, s0_perc = 0.05,
                    group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                    alpha = 0.06, n_cl = 2, verbose = TRUE,
                    global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = F,
                    iter_control = list(tol = 1e-2, max_iter = 100,
                                        verbose = TRUE),
                    em_control = list(tol = 1e-5, max_iter = 100),
                    lme_control = lme4::lmerControl(),
                    mdfdr_control = list(fwer_ctrl_method = Padj_Method, B = 100),
                    trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                nrow = 2,
                                                                byrow = TRUE),
                                                         matrix(c(-1, 0, 1, -1),
                                                                nrow = 2,
                                                                byrow = TRUE)),
                                         node = list(2, 2),
                                         solver = "ECOS",
                                         B = 100))
  
  res = output$res
  
  
  data_ancombc = res %>%
    dplyr::select(taxon, contains(Grup_variable))
  
  ##################
  ##Create plot 1 ###
  ##################


  data_ancombc1 <- data_ancombc[,c(1:7)]
  
  ####Standarize names
  
  colnames(data_ancombc) <- c("Taxon","log2FoldChange","SE","W","p_value","p_adj","Diff")
  

  data_ancombc <-  data_ancombc %>%
    dplyr::mutate(diffexpressed   = case_when(p_value >= 0.05 ~ "No",
                                              p_value <= 0.05 ~ "Yes"
                                              ),
                  delabel = case_when(p_adj >= 0.05 ~ NA,
                                      p_adj <= 0.05 ~ Taxon
                                      )
                  )
  
  
  p <- ggplot(data=data_ancombc, aes(x=log2FoldChange, y=-log10(p_value), col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept= 0, col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
  
  
  ##################
  ##Save the data###
  ##################
  
  
  # write.csv(res,file = paste0(Results_directory,"/",Analysis_name,".csv"))
  # write.xlsx(res, paste0(Results_directory,"/",Analysis_name,".xlsx"))
  # ggsave(plot = p,paste0(Results_directory,"/",Analysis_name,".png"), width = 20, height = 20, units = "cm")
  
  return(list(res,p))
  
}





###################################################
########RUN ANCOMBC(For non source users)##########
###################################################

# Diferential_function(Physeq,Taxonomy,Padj_Method,Adjformula,Grup_variable)


