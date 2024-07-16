###################################################################
###################################################################
###################################################################
###     VALENCA R adapted                                    ######
###  Linkedin: https://www.linkedin.com/in/adria-cruells/    ######
###  ORCID: https://orcid.org/0000-0002-1179-7997            ######
###################################################################
###################################################################



yue_distance <- function(row, median) {
  taxon_count <- 1
  median_times_obs <- numeric(length(row))
  median_minus_obs_sq <- numeric(length(row))
  for (taxon_abund in row) {
    median_times_obs[taxon_count] <- as.numeric(median[taxon_count]) * taxon_abund
    median_minus_obs_sq[taxon_count] <- (as.numeric(median[taxon_count]) - taxon_abund)^2
    taxon_count <- taxon_count + 1
  }
  product <- sum(median_times_obs, na.rm = TRUE)
  diff_sq <- sum(median_minus_obs_sq, na.rm = TRUE)
  yue_med_dist <- product / (diff_sq + product)
  return(yue_med_dist)
}

CSTs <- c('I-A', 'I-B', 'II', 'III-A', 'III-B', 'IV-A', 'IV-B', 'IV-C0', 'IV-C1', 'IV-C2', 'IV-C3', 'IV-C4', 'V')

reference_centroids<- read.csv("./VALENCIA-master/CST_centroids_012920.csv", sep = ',')



sample_data_OG <- readRDS(file =   "./../vaginal.rds")
sample_data_OG <- as.data.frame(otu_table(sample_data_OG))
sample_data_OG$Taxa <- rownames(sample_data_OG)


library(dplyr)

Genus_to_filter <- gsub(pattern = "g_",replacement = "",x =  colnames(reference_centroids))
Genus_to_filter <- gsub(pattern = "f_",replacement = "",x =  Genus_to_filter)

colnames(reference_centroids) <- Genus_to_filter
Genus_to_filter <- Genus_to_filter[-1]


CST_table <- data.frame()

for (i in Genus_to_filter) {
  
  if(length(strsplit(i,split = "_")[[1]]) == 2){
    
    a  <- sample_data_OG %>%
      dplyr::filter(grepl(strsplit(i,split = "_")[[1]][2],Taxa))%>%
      select(-ncol(.))
    
    sum_columnas <- cbind.data.frame(i,t(data.frame(colSums(a))))
    
    CST_table <- rbind.data.frame(CST_table,sum_columnas)
    
  }else{
    a  <- sample_data_OG %>%
      dplyr::filter(grepl(i,Taxa))%>%
      select(-ncol(.))
    
    sum_columnas <- cbind.data.frame(i,t(data.frame(colSums(a))))
    
    CST_table <- rbind.data.frame(CST_table,sum_columnas)
  }
}


colnames(CST_table)[1] <- "TAXA"
rownames(CST_table) <- NULL

CST_table[,-1] <- CST_table[,-1]/colSums(CST_table[,-1])


CST_table <- as.data.frame(t(CST_table))

colnames(CST_table) <- CST_table[1,]
CST_table <- CST_table[-1,]
CST_table$sub_CST <- rownames(CST_table)


combined_data <- rbind(CST_table, reference_centroids)
sample_data <- combined_data[1:(nrow(combined_data) - 13), ]
sample_data[is.na(sample_data)] <- 0
sample_data <- sample_data[, !names(sample_data) %in% 'sub_CST']
reference_centroids <- combined_data[(nrow(combined_data) - 12):nrow(combined_data), ]
reference_centroids <- reference_centroids[, !names(reference_centroids) %in% c("sampleID", "read_count")]
rownames(reference_centroids) <- reference_centroids$sub_CST

sample_data_rel <- sample_data
sample_data_rel <- as.data.frame(apply(sample_data_rel, 2, as.numeric))



Results <- data.frame(c(rownames(sample_data)))

for (CST in CSTs) {
  Results[paste0(CST, '_sim')] <- apply(sample_data_rel, 1, function(x) yue_distance(x, reference_centroids[CST, ]))
}




columna_con_max <- function(row) {
  nombres_columnas <- colnames(Results)
  indice_max <- which.max(row)
  nombre_columna_max <- nombres_columnas[indice_max]
  return(nombre_columna_max)
}


nombres_columnas_max <- apply(Results, 1, columna_con_max)
Results_Final <- cbind.data.frame(Results,nombres_columnas_max)

colnames(Results_Final) <- c("SampleID","I-A_sim","I-B_sim","II_sim","III-A_sim","III-B_sim",
                              "IV-A_sim","IV-B_sim","IV-C0_sim","IV-C1_sim","IV-C2_sim","IV-C3_sim",               
                              "IV-C4_sim","V_sim","CST_group" )


write.csv(Results_Final, paste0('CST_Valencia.csv'), row.names = FALSE)

sample_data_OG2 <- readRDS(file =   "./../vaginal.rds")

sample_data_OG2 <- data.frame(sample_data(sample_data_OG2))

colnames(sample_data_OG2)[2] <- "SampleID"


Results_Meta <- merge(Results_Final,sample_data_OG2,by = "SampleID")

write.csv(Results_Final, paste0('CST_Valencia_META.csv'), row.names = FALSE)



library(plyr)

Graph_table <- plyr::rbind.fill(data.frame(table(Results_Meta$CST_group[Results_Meta$Pacient == "Yes"])/sum(table(Results_Meta$CST_group[Results_Meta$Pacient == "Yes"]))),
                                data.frame(table(Results_Meta$CST_group[Results_Meta$Pacient == "No"])/sum(table(Results_Meta$CST_group[Results_Meta$Pacient == "No"]))))

colnames(Graph_table) <- CSTs[-10] ### No hay C2

Graph_table$Pacient <- c("Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes","Yes",
                         "Yes","No","No","No","No","No","No","No","No","No","No")



library(ggplot2)


ggplot(Graph_table,aes(x = Pacient,y = Freq,fill = Var1)) + geom_col()




####Statistical analysis
library(tidyr) 
library(MASS) 


Graph_table2 <- Graph_table %>%
  pivot_wider(names_from = Pacient, values_from = Freq)




colnames(Graph_table2) <- c("CST","Pacient","Control")



library(freqtables)
library(dplyr)

test_table <- Results_Meta %>%
  dplyr::select(CST_group,Pacient)%>%
  freq_table(CST_group,Pacient)

write.csv(test_table, paste0('CST_Valencia_META.csv'), row.names = FALSE)

#### perfrom the analysis by the sexual effect



Results_Final <- read.csv("CST_Valencia.csv")

Meta <- data.frame(sample_data(Vaginal))


Results_Final <- merge(Results_Final,Meta,by.x = "SampleID",by = "colnames.mostres_filt.")

test_table <- Results_Final %>%
  mutate(
    Post_coital = case_when(Pacient == "No" ~ "Control",
                            Pacient == "Yes" &  Postcoital.UTI == "Yes"~ "UTI-postcoital",
                            Pacient == "Yes" &  Postcoital.UTI == "No"  ~ "UTI-no-postcoital")
  )%>%
  dplyr::select(Post_coital,CST_group)%>%
  freq_table(Post_coital,CST_group) %>%
  mutate(
    Lacto_dominated = case_when(col_cat == "I-A_sim" ~ "Lactobacillus dominated",
                            col_cat == "I-B_sim" ~ "Lactobacillus dominated",
                            col_cat == "II_sim" ~ "Lactobacillus dominated",
                            col_cat == "III-A_sim" ~ "Lactobacillus dominated",
                            col_cat == "III-B_sim" ~ "Lactobacillus dominated",
                            col_cat == "IV-A_sim" ~ "Non Lactobacillus dominated",
                            col_cat == "IV-B_sim" ~ "Non Lactobacillus dominated",
                            col_cat == "IV-C0_sim" ~ "Non Lactobacillus dominated",
                            col_cat == "IV-C1_sim" ~ "Non Lactobacillus dominated",
                            col_cat == "IV-C3_sim" ~ "Non Lactobacillus dominated",
                            col_cat == "IV-C4_sim" ~ "Non Lactobacillus dominated",
                            col_cat == "V_sim" ~ "Lactobacillus dominated"))%>%
  mutate(
    Predominant_group = case_when(col_cat == "I-A_sim" ~ "L.crispatus",
                                col_cat == "I-B_sim" ~ "L.crispatus",
                                col_cat == "II_sim" ~ "L.gasseri",
                                col_cat == "III-A_sim" ~ "L. iners",
                                col_cat == "III-B_sim" ~ "L. iners",
                                col_cat == "IV-A_sim" ~ "Candidatus Lachnocurva vaginae",
                                col_cat == "IV-B_sim" ~ "G. vaginalis",
                                col_cat == "IV-C0_sim" ~ "Prevotella genus",
                                col_cat == "IV-C1_sim" ~ "Streptococcus genus",
                                col_cat == "IV-C3_sim" ~ "Bifidobacterium dominated",
                                col_cat == "IV-C4_sim" ~ "Staphylococcus dominated",
                                col_cat == "V_sim" ~ "L. jensenii"))
    


Fig4_F <- ggplot(test_table,aes(row_cat,percent_row, fill = Predominant_group)) + geom_col() + facet_wrap(.~Lacto_dominated) +
  theme(legend.title = element_blank(),legend.position = "bottom")+ scale_fill_brewer(palette = "Paired") + ylab("") + xlab("")+
  theme_bw()


