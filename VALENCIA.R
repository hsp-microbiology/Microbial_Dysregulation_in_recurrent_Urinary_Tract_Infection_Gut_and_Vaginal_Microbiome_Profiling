###################################################################
###################################################################
###################################################################
###     VALENCA R adapted                                    ######
###  Linkedin: https://www.linkedin.com/in/adria-cruells/    ######
###  ORCID: https://orcid.org/0000-0002-1179-7997            ######
###################################################################
###################################################################




#VALENCIA
###El propósito de esta herramienta es clasificar las comunidades microbianas vaginales en tipos de estado comunitario (CST)
###de manera estandarizada y repetible, los grupos se basan en los datos de más de 13,000 mujeres
###Se consideran un total de 13 CSTs: CST I-A, I-B, II, III-A, III-B, IV-A, IV-B, IV-C0, IV-C1, IV-C2, IV-C3, IV-C4, V 
####Este script prueba nuevas muestras basadas en su similitud con los centroides de los grupos definidos previamente

#definiendo la función para determinar theta de Yue-Clayton
yue_distance <- function(row, median) {
  #creando una variable de conteo para indexar la lista mediana
  taxon_count <- 1
  #creando listas para almacenar los resultados iterativamente
  median_times_obs <- numeric(length(row))
  median_minus_obs_sq <- numeric(length(row))
  #recorriendo la fila y calculando el producto y la diferencia al cuadrado entre los datos de la fila y los datos de la mediana
  for (taxon_abund in row) {
    #calculando p * q
    median_times_obs[taxon_count] <- as.numeric(median[taxon_count]) * taxon_abund
    #calculando (p-q) al cuadrado
    median_minus_obs_sq[taxon_count] <- (as.numeric(median[taxon_count]) - taxon_abund)^2
    taxon_count <- taxon_count + 1
  }
  #calculando la suma de p * q
  product <- sum(median_times_obs, na.rm = TRUE)
  #calculando la suma de (p-q) al cuadrado
  diff_sq <- sum(median_minus_obs_sq, na.rm = TRUE)
  #calculando la distancia de Yue
  yue_med_dist <- product / (diff_sq + product)
  #retornando el valor de la distancia de Yue
  return(yue_med_dist)
}

#lista de subCSTs 
CSTs <- c('I-A', 'I-B', 'II', 'III-A', 'III-B', 'IV-A', 'IV-B', 'IV-C0', 'IV-C1', 'IV-C2', 'IV-C3', 'IV-C4', 'V')

#leyendo los centroides de los CSTs de entrada
reference_centroids<- read.csv("./VALENCIA-master/CST_centroids_012920.csv", sep = ',')
#leyendo la tabla de muestras que se probarán contra los centroides


sample_data_OG <- readRDS(file =   "./../vaginal.rds")
sample_data_OG <- as.data.frame(otu_table(sample_data_OG))
sample_data_OG$Taxa <- rownames(sample_data_OG)
###Get the only taxas used to obtain the Vaginal types

library(dplyr)
## Obtain all genus

Genus_to_filter <- gsub(pattern = "g_",replacement = "",x =  colnames(reference_centroids))
Genus_to_filter <- gsub(pattern = "f_",replacement = "",x =  Genus_to_filter)

colnames(reference_centroids) <- Genus_to_filter
Genus_to_filter <- Genus_to_filter[-1]

### Sum all the bacteria included in each superior level, such as Families

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

###Make relatives

CST_table[,-1] <- CST_table[,-1]/colSums(CST_table[,-1])



### PRepare the table

CST_table <- as.data.frame(t(CST_table))

colnames(CST_table) <- CST_table[1,]
CST_table <- CST_table[-1,]
CST_table$sub_CST <- rownames(CST_table)





#forzando que los datos de muestra tengan el mismo número de columnas que la referencia y en el mismo orden
combined_data <- rbind(CST_table, reference_centroids)
sample_data <- combined_data[1:(nrow(combined_data) - 13), ]
sample_data[is.na(sample_data)] <- 0
sample_data <- sample_data[, !names(sample_data) %in% 'sub_CST']
reference_centroids <- combined_data[(nrow(combined_data) - 12):nrow(combined_data), ]
reference_centroids <- reference_centroids[, !names(reference_centroids) %in% c("sampleID", "read_count")]
rownames(reference_centroids) <- reference_centroids$sub_CST

#convirtiendo todas las cuentas de lectura a datos de abundancia relativa y agregando las dos primeras columnas nuevamente
sample_data_rel <- sample_data
sample_data_rel <- as.data.frame(apply(sample_data_rel, 2, as.numeric))


#midiendo la similitud de cada muestra con cada centroide de subCST utilizando theta de Yue + Clayton

Results <- data.frame(c(rownames(sample_data)))

for (CST in CSTs) {
  Results[paste0(CST, '_sim')] <- apply(sample_data_rel, 1, function(x) yue_distance(x, reference_centroids[CST, ]))
}

###Check the most accurated 


columna_con_max <- function(row) {
  nombres_columnas <- colnames(Results)
  indice_max <- which.max(row)
  nombre_columna_max <- nombres_columnas[indice_max]
  return(nombre_columna_max)
}

# Aplicar la función a lo largo de las filas del data frame
nombres_columnas_max <- apply(Results, 1, columna_con_max)
Results_Final <- cbind.data.frame(Results,nombres_columnas_max)

colnames(Results_Final) <- c("SampleID","I-A_sim","I-B_sim","II_sim","III-A_sim","III-B_sim",
                              "IV-A_sim","IV-B_sim","IV-C0_sim","IV-C1_sim","IV-C2_sim","IV-C3_sim",               
                              "IV-C4_sim","V_sim","CST_group" )

#exportar las asignaciones en un nuevo archivo CSV
write.csv(Results_Final, paste0('CST_Valencia.csv'), row.names = FALSE)

sample_data_OG2 <- readRDS(file =   "./../vaginal.rds")

sample_data_OG2 <- data.frame(sample_data(sample_data_OG2))

colnames(sample_data_OG2)[2] <- "SampleID"


Results_Meta <- merge(Results_Final,sample_data_OG2,by = "SampleID")

write.csv(Results_Final, paste0('CST_Valencia_META.csv'), row.names = FALSE)


###Make an example graph
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


