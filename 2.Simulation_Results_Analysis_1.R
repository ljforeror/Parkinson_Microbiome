############################################Libraries used


library("tidyverse")
library("ggplot2")
library("data.table")
library("sjPlot")
library("rstatix")
library("reshape2")
library("ggpubr")
library("ggplot2")
library("qqplotr")
library("parallel")
library("coin")
library("viridis")
library("RcmdrMisc")
library("plot.matrix")
library("patchwork")
library("nlme")
library("Matrix")
library("lmerTest")
library("lme4")
library("kableExtra")
library("hrbrthemes")
library("Hmisc")
library("gtable")
library("gridGraphics")
library("gridExtra")
library("GGally")
library("egg")
library("dplyr")
library("datasets")
library("cowplot")
library("corrplot")
library("broom.mixed")
library("amap") #para Dist
library("vtable")
library("vegan")
library("tidyr")
library("survival")
library("stargazer")
library("reshape")
library("readxl")
library("readr")
library("phyloseq")
library("metagMisc")
library("gridExtra")
library("grid")
library("ggplotify")
library("fantaxtic")
library("factoextra")
library("doBy")
library("DESeq2")
library("dendextend")
library("car")
library("ape")
library(ghibli)
library("DescTools")

####################################### Data bases


metadata <- read.table("F:/articulo/decompri/Metadata_english.txt")

# Computational modelling analysis

##Pattern list 
patron.list1 <- vector("list",4)
patron.list1[[1]] <- "concentraciones.tsv"
patron.list1[[2]] <- "flujos.tsv"
patron.list1[[3]] <- "Biomasa.tsv"
patron.list1[[4]] <- "Celulas.tsv"

patron.list2 <- vector("list",3)
patron.list2[[1]] <- "flujosEspecieHora.tsv"
patron.list2[[2]] <- "flujos2EspecieHora.tsv"
patron.list2[[3]] <- "flujosTotalesPorEspeciePorHora.tsv"


#################################################################################################
###################################################################################################
#################################################################################################

###############1. Concentrations

for (n in 1:2) { 
  dir.rds <- dir(path = "F:/articulo/segundo_1",pattern = "concentraciones.tsv" ,full.names = TRUE)
  
  #get a list with each pattern from each individuo
  lista1 <- vector("double") 
  for(individuo in 1:length(dir.rds)){  
    h1 <- read.table(file = dir.rds[[individuo]], header = TRUE)
    h1$nombre <- dir.rds[[individuo]]
    #get the list of file from each pattern
    lista1[[individuo]] <- list(h1)
    names(lista1)[[individuo]] <- dir.rds[[individuo]]
  }
  
  #merge the list to get one from all the individuo
  #https://blog.zhaw.ch/datascience/r-reduce-applys-lesser-known-brother/
  #merge in a way that you do not lose info (all=T)
  patron.list1[[1]] <- Reduce(function(d1,d2)merge.data.frame(d1,d2,all = T), lista1)
}
concentraciones_junio <- as.data.frame(patron.list1[[1]])

concen<-concentraciones_junio %>% 
  mutate_if(is.character, str_replace, 'F:/articulo/segundo_1/Healthy_', '')

concen1<-concen %>% 
  mutate_if(is.character, str_replace, 'F:/articulo/segundo_1/Patient_', '')

concen2<-concen1%>% 
  mutate_if(is.character, str_replace, '_.RDS-concentraciones.tsv', '')


####
aqw <- reshape2::dcast(concen2, time + replc + nombre ~ sub, value.var= "value")
######## correction by biomass
concentration_df <- aqw

#Changing NAs by cero
concentration_df[is.na(concentration_df)] <- 0
#taking tables for each hour of interest
concentration_df_new <- concentration_df
concentration_df_1b <- concentration_df[concentration_df["time"]==2,]
concentration_df_12b <- concentration_df[concentration_df["time"]==12,]

#biomass hour 12 less hour 2
concentration_df_hb <- concentration_df_12b[,"biomass"]-concentration_df_1b[,"biomass"]

#correcting metabolites 12h concentration by the biomass
for (i in 4:ncol(concentration_df_12b)){
    concentration_df_12b[,i] <- (concentration_df_12b[,i] * concentration_df_hb /concentration_df_12b[,"biomass"])
                     }

######### Correcting by time
#table for each time
concentration_df_1 <- concentration_df[concentration_df["time"]==1,]
concentration_df_12 <- concentration_df_12b[concentration_df_12b["time"]==12,]

#metabolites concentration 12h less 1h
concentration_df_h <- concentration_df_12[,-c(1:3)]-concentration_df_1[,-c(1:3)]

#eliminating variables without importance and organizing table for analysis
ncol(concentration_df_h)
concentration_df_h$numero <- c(1:500)
concentration_df_12$numero <- c(1:500)
uji <- full_join(concentration_df_12[,c(1:3, 139)],abs(concentration_df_h))

#Erase column
concentration_df_12 <- uji[,-c(4)]

###################################Take off metabolites with too many cero values (ceros)
porcentaje <- colSums(concentration_df_12[, 4:ncol(concentration_df_12)]==0)/nrow(concentration_df_12)
#percentage > 0.9
concentration_df_12_f <- concentration_df_12[,c(rep(TRUE,3), !(porcentaje > 0.9))]


#############################reshape table

r <- reshape2::melt(concentration_df_12, id.var = c("time", "replc", "nombre"))

###########################mean replicates
b <- aggregate(value~variable+time+nombre,r,mean)
colnames(r)
#########################bind with  metadata
###column names columns
colnames(metadata)[1] <- "nombre"
##changing some data

metadata$Sex <- gsub( "M", "F", metadata$Sex)
metadata$Sex <- gsub( "H", "M", metadata$Sex)

metadata$Alcohol <- gsub( "SI", "Yes", metadata$Alcohol)
metadata$Smoke <- gsub( "SI", "Yes", metadata$Smoke)

metadata$Disease <- metadata$nombre
metadata$nombre <- metadata$ID

metadata$nombre <- gsub( "c", "H", metadata$nombre)
metadata$nombre <- gsub( "p", "P", metadata$nombre)
metadata$ID <- NULL

metadata$Parkinson <- gsub( "0", "Healthy", metadata$Parkinson)
metadata$Parkinson <- gsub( "1", "Parkinson", metadata$Parkinson)

colnames(metadata)[colnames(metadata) == 'PE'] <- 'PD'

#taked/modified from: https://www.iteramos.com/pregunta/73573/en-la-r-como-reemplazo-el-texto-dentro-de-una-cadena
Tabla_real_10<-full_join(b,metadata)
Tabla_real_10<-Tabla_real_10[Tabla_real_10$variable!="biomass", ]

######guardar
########################reshape
#Tabla_real_10[["scaled_2"]] <- log(Tabla_real_10$value)
Tabla_real_10_Ancha <- reshape2::dcast(Tabla_real_10, 
                                       nombre ~ variable,
                                       value.var = "value")

##########################################################12 h GLM ##################################
########################get the table for each metabolite

li <- as.character(unique(Tabla_real_10$variable))
li
list <- vector(mode = "list",length = length(li))

for ( j in 1:length(li) ){  
  
  d <- Tabla_real_10[which(Tabla_real_10$variable==li[j]),]#c (1:11,24,35,43:46)]
  
  list[[j]] <-d
  names(list)[[j]] <- li[j]
  
}


#####
################get  each of the tables and filter by the 12 hours
list2 <- vector(mode = "list",length = length(list))

for ( e in 1:length(list) ){  
  
  h1 <- as.data.frame(list[[e]])
  list2[[e]] <-h1
  
  names(list2)[[e]] <- names (list)[[e]]
  
  
}


###################for each table with each metabolite (12H) get the glmer


list3 <- vector(mode = "list",length = length(list2))

for ( i in 1:length(list2) ){  
  #i=2
  h1 <- as.data.frame(list2[[i]])
  h1$Parkinson <- as.factor(h1$Parkinson)
  
  if (length (unique(h1$nombre)) == 1) {
    list3[[i]]=NULL
  }
  else {
    
    possibleError <- tryCatch({
      
      p <- glm(data=h1, factor(Parkinson) ~ value, family = "binomial")
      
      
      list3[[i]] <-p
    },
    error=function(e) {
      e
      return(print("no se guarda este modelo"))}
    )
    
    if(inherits(possibleError, "error")){
      
      #print("ERROR:)")
      list3[[i]]=NULL
    }
    
  }
  
  names(list3)[[i]] <- names (list2)[[i]]
}


#######################save the coeeficients in a list
list4 <- vector(mode = "list",length = length(list3))

for ( t in 1:length(list3) ){  
  
  
  if (length (unique(h1$nombre)) == 1) {
    list4[[t]]=NULL
  }
  else {
    
    possibleError <- tryCatch( {
      
      x <- coef(summary(list3[[t]]))
      list4[[t]] <-x
    },
    error=function(e) {
      e
      return(print("no se guarda este modelo"))}
    
    )
    
    if(inherits(possibleError, "error")){
      
      #print("ERROR:)")
      list4[[t]]=NULL
    }
    
  }
  
  ####
  
  names(list4)[[t]] <- names (list3)[[t]]
}




###############################get the second coeficient and the merge them in a dataframe

patron <- vector(mode = "list",length = 3)

#get a list with each pattern from each individuo
lista7 <- vector("double") 
for (n in 1:2) { 
  for(individuo in 1:length(list4)){ 
    
    h1 <- as.data.frame(list4[[individuo]], header = TRUE)
    if (length(h1)== 0) {
      
    }
    
    h1$ID <-names(list4)[[individuo]]##
    #get the list of file from each pattern
    lista7[[individuo]] <- list(h1[2,])
    
    names(lista7)[[individuo]] <- names (list4)[[individuo]]
  } 
  patron[[1]] <- Reduce(function(d1,d2)merge.data.frame(d1,d2,all = T), lista7)
  
}


rbind(unlist(lista7[[1]]), unlist(lista7[[133]]))

######## merge al the list in a data.frame

o <- do.call("rbind", lista7)
q2 <- Reduce(function(d1,d2)merge.data.frame(d1,d2,all = T), o)

######################################################correction of p-values 

q2$FDR <- p.adjust(q2[,4],method="fdr")
q2$bon <- p.adjust(q2[,4],method="bonferroni", n = length(q2[,4]))
q2$holm <- p.adjust(q2[,4],method="holm")
q2$hochberg <- p.adjust(q2[,4],method="hochberg")
q2$hommel <- p.adjust(q2[,4],method="hommel")
q2$BH <- p.adjust(q2[,4],method="BH")
q2$BY <- p.adjust(q2[,4],method="BY")
q2$none <- p.adjust(q2[,4],method="none")


#check in a different table with p.value< 0.05 
solo_q2<- q2[c(which(q2[,4]<=0.05)),]##### pvalue non corrected
list_q2<- list(solo_q2$ID)
#saving info
write.table(format(solo_q2, scientific=FALSE, digits=4),file ="C:/Users/57322/Desktop/cod/gml_concentration12H.csv", dec=",")




## PCA

################################PCA metabolites
Healthy <-  c(rep("Healthy",25))
Patient <- c(rep("Patient",25))
rownames(Tabla_real_10_Ancha) <- Tabla_real_10_Ancha$nombre
Tabla_real_10_Ancha_1 <- Tabla_real_10_Ancha[,-1]
Tabla_real_10_Ancha$Parkinson <- c(Healthy, Patient)
Tabla_real_10_Ancha_1 <- Tabla_real_10_Ancha[,2:134]

res.pca_1 <- prcomp(Tabla_real_10_Ancha_1, scale = TRUE)
fviz_eig(res.pca_1)


pca_met2 <-fviz_contrib(res.pca_1, choice = "var", axes = 2, top = 25)
pca_met3 <-fviz_contrib(res.pca_1, choice = "var", axes = 1, top = 25)

fviz_pca_biplot(res.pca_1, axes = c(1, 2), addEllipses = TRUE, repel = TRUE, select.var = list(contrib = 20),  label="var",
                col.var = "#63B8FF", # Variables color
                col.ind = groups,  # Individuals color
                palette = c("#CD6889", "#8B668B"))

fviz_screeplot(res.pca_1, addlabels = TRUE, ylim = c(0, 20))

groups <- as.factor(Tabla_real_10_Ancha$Parkinson)

pca_met <- fviz_pca_ind(res.pca_1,
                        col.ind = groups, # color by groups
                        palette = c("#CD6889", "#8B668B"),
                        addEllipses = TRUE, # Concentration ellipses
                        ellipse.type = "confidence",
                        legend.title = "Groups",
                        repel = TRUE
)
#saving info

#ggsave(pca_met,                                        #nombre de la gr??fica en R
#file=paste(resultado,"pca_met_metabolites", ".jpeg", sep=''), #Nombre del jpeg
#height=6, width=5, units="in", dpi=300, bg = "white")
#ggsave(pca_met1,                                        #nombre de la gr??fica en R
#file=paste(resultado,"pca_met1_metabolites", ".jpeg", sep=''), #Nombre del jpeg
#height=6, width=5, units="in", dpi=300, bg = "white")
#ggsave(pca_met2,                                        #nombre de la gr??fica en R
#file=paste(resultado,"pca_met2_metabolites", ".jpeg", sep=''), #Nombre del jpeg
#height=6, width=5, units="in", dpi=300, bg = "white")
#ggsave(pca_met3,                                        #nombre de la gr??fica en R
#file=paste(resultado,"pca_met3_metabolites", ".jpeg", sep=''), #Nombre del jpeg
#height=6, width=5, units="in", dpi=300, bg = "white")





##############################spearmann correlation

listcor <- vector(mode = "list",length = length(list))
variablesitas <- c("value", "Peso", "Fuma", "Alcohol","deposicionessemana")

for ( t in 1:length(list)){  
  #print(t)
  h3 <- as.data.frame(list[[t]])
  for ( var in 1:length(variablesitas)){
    cores <- cor.test(as.numeric(h3$Scale), as.numeric(h3$value),
                      method="spearman")
   
    
    listcor[[t]] <-c(cores$p.value, cores$estimate)
    
    names(listcor)[[t]] <- names (list)[[t]]
  }
}


#summarie important information
ooq <- do.call("rbind", listcor)
qcor <- as.data.frame(ooq)
qcor$FDR <- p.adjust(qcor[,1],method="fdr")
qcor[abs(qcor$rho) > 0.4,]





#bars and violin figures Important metabolites and scfas(short chain fatty acids)

#prepare table for figures

concentration_df_12_4 <- concentration_df_12_f


Healthy <-  c(rep("Healthy",25))
Patient <- c(rep("Patient",25))
df <- concentration_df_12_f[1:50,]

df$PD <- c(Healthy, Patient)
r4 <- reshape2::melt(df, id.var = c("time", "replc", "nombre","PD" ))
b4 <- aggregate(value~time+nombre+PD+variable,r4,mean)
b5 <- b4
#prepare variable names
new_variables <- c("Phenylacetic acid", "Indole","L_Tryptophan", "Fructose", "Myristic acid", "N-Acetylneuraminic acid", "3-Methyl-2-Oxovaleric Acid")
ex_reac <- c("pac", "indole", "trp_L", "fru", "ttdca", "acnam","3mop")
shorts<- c("ac", "but", "ppa")
Shorts_new <- c("Acetate" , "Butyrate", "Propionate")

#taken from: https://stackoverflow.com/questions/19424709/r-gsub-pattern-vector-and-replacement-vector
for (i in seq_along(ex_reac)){
  b4$variable <- gsub(ex_reac[i], new_variables[i], b4$variable)
}

#### 7 important metabolites
listaazucar <-  vector("list")#, length = 10)
listviolin <-  vector("list")
for(i in new_variables){
  #i = "acetate"
  #print(i)
  
  
  hc <- b4[which(b4$variable == i),]
  
  df.mean = hc %>% 
    group_by(PD) %>% 
    mutate(ymean = mean(value))
  
  bars_metabol_plot <-ggplot(hc, aes(x = nombre, y = value, fill = PD)) +
    geom_col(width = 0.9,position=position_dodge(0.9)) +
    labs(title=paste(i)) +
    coord_flip() +
    geom_errorbar(data=df.mean, aes(nombre, ymax = ymean, ymin = ymean),
                  size=1, linetype = "longdash", inherit.aes = F, width = 1, color = "cadetblue1")+
    scale_fill_manual(values = c("#CD6889", "#8B668B"))+
    scale_y_continuous(limits = c(0,NA), expand = c(0, 0)) +
    labs( 
      x = "",
      y = "Concentration mM", xlab(i)) +
   
  theme_classic()
  graph<-bars_metabol_plot +
    theme(axis.text.x = element_text(face="plain", color="black", 
                                                     size=16, angle=0),
                          axis.text.y = element_text(face="plain", color="black", 
                                                     size=16, angle=0, hjust = 1), 
                          axis.title.x = element_text(face="plain", color="black",
                                                    size=16, angle=0, hjust = 1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'black'), legend.position="bottom", axis.title.y = element_blank(), axis.title.x =  element_text(hjust = 0.5))
 
  
  plot(graph)
  ggsave(graph,                                        #name 
         file=paste("C:/Users/57322/Desktop/cod/",i , ".jpeg", sep='_'), #Name jpeg
         height = 12, width = 14, dpi=300, bg = "white") 
  

  listaazucar[[i]] <- graph

  ####violin
  violin_graphs <-ggplot(data = hc,aes(x = variable, y = value, fill = PD))+
    #scale_fill_viridis_d( option = "D")+  
    geom_violin(alpha=0.7, position = position_dodge(width = .75),size=1,color="black") +      
    #geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+      
    geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
    scale_fill_manual(values = c("#CD6889", "#8B668B"))+
    labs( 
      x = "",
      y = "Concentration mM", xlab(i)) +
    theme(axis.text.x = element_text(angle = 0)) +
  theme_pubr()
  violin_graphs<-violin_graphs + theme(axis.text.x = element_text(face="plain", color="black", 
                                                                       size=16, angle=0),
                                            axis.text.y = element_text(face="plain", color="black", 
                                                                       size=16, angle=0)) +
    theme(legend.position="bottom" )
  
  
  plot(violin_graphs)
  ggsave(violin_graphs,                                        
         file=paste("C:/Users/57322/Desktop/cod/",i , ".jpeg", sep='_'), 
         height = 12, width = 14, dpi=300, bg = "white") 
  
  
  listviolin[[i]] <- violin_graphs
  
}

yleft_plo <- textGrob(expression(paste("Concentration mM")), 
                  rot = 90, gp = gpar(fontsize = 20))
xleft_plo <- textGrob(expression(paste("Concentration mM")), 
                      rot = 0, gp = gpar(fontsize = 20), vjust = -4)



combined =  listaazucar[[1]] + listaazucar[[2]] + listaazucar[[3]] + listaazucar[[4]] + listaazucar[[5]] + listaazucar[[6]] + listaazucar[[7]] & theme(legend.position = "bottom" , axis.title.x = element_blank())
Plot_metab_10<- combined + plot_layout(guides = "collect", ncol = 7)
dbars<-grid.arrange(patchworkGrob(Plot_metab_10), bottom= xleft_plo)
#saving info
ggsave("C:/Users/57322/Desktop/cod/barras_metabol.png", plot = dbars, width = 25, height = 15, units = "in", dpi = 300)

combined =  listviolin[[1]]  + listviolin[[2]] + listviolin[[3]] + listviolin[[4]] + listviolin[[5]] + listviolin[[6]] + listviolin[[7]] & theme(legend.position = "bottom", axis.title.y = element_blank())
Plot_metab_10<-  combined + plot_layout(guides = "collect", ncol = 7)
dd<-grid.arrange(patchworkGrob(Plot_metab_10), left = yleft_plo)
#saving info
ggsave("C:/Users/57322/Desktop/cod/violin_metabol.png", plot = dd, width = 25, height = 15, units = "in", dpi = 300)

##################short chain fatty acids

for (i in seq_along(shorts)){
  b4$variable <- gsub(shorts[i], Shorts_new[i], b4$variable)
} 
###vectors list
listaazucars <-  vector("list")#, length = 10)
listviolins <-  vector("list")
for(i in Shorts_new){
  
  
  
  hc <- b4[which(b4$variable == i),]
  
  df.mean = hc %>% 
    group_by(PD) %>% 
    mutate(ymean = mean(value))
  
  bars_metabol_plot <-ggplot(hc, aes(x = nombre, y = value, fill = PD)) +
    geom_col(width = 0.9,position=position_dodge(0.9)) +
    labs(title=paste(i)) +
    coord_flip() +
    geom_errorbar(data=df.mean, aes(nombre, ymax = ymean, ymin = ymean),
                  size=1, linetype = "longdash", inherit.aes = F, width = 1, color = "cadetblue1")+
    scale_fill_manual(values = c("#CD6889", "#8B668B"))+
    scale_y_continuous(limits = c(0,NA), expand = c(0, 0)) +
    labs( 
      x = "",
      y = "Concentration mM", xlab(i)) +
    
    theme_classic()
  graph<-bars_metabol_plot +
    theme(axis.text.x = element_text(face="plain", color="black", 
                                     size=14, angle=0),
          axis.text.y = element_text(face="plain", color="black", 
                                     size=14, angle=0, hjust = 1), 
          axis.title.x = element_text(face="plain", color="black",
                                      size=14, angle=0, hjust = 1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'black'), legend.position="bottom", axis.title.y = element_blank(), axis.title.x =  element_text(hjust = 0.5))
  
  
  plot(graph)
  ggsave(graph,                                        
         file=paste("C:/Users/57322/Desktop/cod/",i , ".jpeg", sep='_'), 
         height = 12, width = 14, dpi=300, bg = "white") 
  
  
  listaazucars[[i]] <- graph

  ####violin
  violin_graphs <-ggplot(data = hc,aes(x = variable, y = value, fill = PD))+
     
    geom_violin(alpha=0.7, position = position_dodge(width = .75),size=1,color="black") +      
      geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
      scale_fill_manual(values = c("#CD6889", "#8B668B"))+
      labs( 
      x = "",
      y = "Concentration mM", xlab(i)) +
    theme(axis.text.x = element_text(angle = 0)) +
    theme_pubr()
  violin_graphs<-violin_graphs + theme(axis.text.x = element_text(face="plain", color="black", 
                                                                  size=14, angle=0),
                                       axis.text.y = element_text(face="plain", color="black", 
                                                                  size=14, angle=0)) +
    theme(legend.position="bottom" )
  
  
  plot(violin_graphs)
  ggsave(violin_graphs,                                        
         file=paste("C:/Users/57322/Desktop/cod/",i , ".jpeg", sep='_'), 
         height = 12, width = 14, dpi=300, bg = "white") 
  
  
  listviolins[[i]] <- violin_graphs
  
}



#taked from :https://patchwork.data-imaginist.com/articles/guides/layout.html
combined =  listaazucars[[1]] + listaazucars[[2]] + listaazucars[[3]] & theme(legend.position = "bottom" , axis.title.x = element_blank())
Plot_metab_10<- combined + plot_layout(guides = "collect", ncol = 3)
dbars<-grid.arrange(patchworkGrob(Plot_metab_10), bottom= xleft_plo)
#save result
ggsave("C:/Users/57322/Desktop/cod/barras_metabol_scafs.png", plot = dbars, width = 25, height = 15, units = "in", dpi = 300)

combined =  listviolins[[1]]  + listviolins[[2]] + listviolins[[3]] & theme(legend.position = "bottom", axis.title.y = element_blank())
Plot_metab_10<-  combined + plot_layout(guides = "collect", ncol = 3)
#taked from : https://stackoverflow.com/questions/65291723/merging-two-y-axes-titles-in-patchwork
dd<-grid.arrange(patchworkGrob(Plot_metab_10), left = yleft_plo)
#Save result
ggsave("C:/Users/57322/Desktop/cod/violin_metabol_scafs.png", plot = dd, width = 25, height = 15, units = "in", dpi = 300)





##########################################fin
