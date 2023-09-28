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
library("ghibli")
library("DescTools")
library("ggrepel")
library(wesanderson)


####################################### Data bases


metadata <- read.table("F:/articulo/decompri/Metadata_english.txt")

##fluxes_especies


#pattern list 
## Cells number

patron.list1 <- vector("list",6)
patron.list1[[1]] <- "concentraciones2Hora.tsv"
patron.list1[[2]] <- "concentracionesHora.tsv"
patron.list1[[3]] <- "concentracionesTotalesHora.tsv"
patron.list1[[4]] <- "Celulas.tsv"
patron.list1[[5]] <- "flujosHora.tsv"
patron.list1[[6]] <- "flujosTotalesHora.tsv"

patron.list2 <- vector("list",3)
patron.list2[[1]] <- "flujosEspecieHora.tsv"
patron.list2[[2]] <- "flujos2EspecieHora.tsv"
patron.list2[[3]] <- "flujosTotalesEspecieHora.tsv"
###
#read the pattern in the folder (inside of the folder there are tables from 25 ptes and 25 controls)
#aim: put all the files with the same pattern together in a list, and merge.


#######
rxn.test.dt <- data.table()
for (n in 1:2) { 
  dir.rds <- dir(path = "F:/articulo/segundo_2/",pattern = "flujosTotalesEspecieHora.tsv" ,full.names = TRUE)
  #get a list with each pattern from each individuo
  lista1 <- vector("double") 
  for(individuo in 1:length(dir.rds)){  
    
    h1 <- read.table(file = dir.rds[[individuo]], header = TRUE)
    h1$Sample <- dir.rds[[individuo]]
    #get the list of file from each pattern
    lista1[[individuo]] <- list(h1)
    
  }
  
  #merge the list to get one from all the individuo
  #https://blog.zhaw.ch/datascience/r-reduce-applys-lesser-known-brother/
  #merge in a way that you do not lose info (all=T)
  patron.list2[[1]] <- Reduce(function(d1,d2)merge.data.frame(d1,d2,all = T), lista1)
  rxn.test.dt <- rbind(rxn.test.dt, h1)
}
species_fluxes_2022 <- as.data.frame(patron.list2[[1]])
species_flujos <- species_fluxes_2022[,c(1:4)]
species_flujos$mmol_per_gram_dry_weight_per_hour  <- rowSums(species_fluxes_2022[,5:54], na.rm = TRUE )


sflux<-species_flujos %>% 
  mutate_if(is.character, str_replace, 'F:/articulo/segundo_2/Healthy_', '')

sflux1<-sflux %>% 
  mutate_if(is.character, str_replace, 'F:/articulo/segundo_2/Patient_', '')


sflux2<-sflux1%>% 
  mutate_if(is.character, str_replace, '_.RDS-flujosTotalesEspecieHora.tsv', '')

#joha_12 <- www[www["time"]==12,]
####
wq<-sflux2[with(sflux2, order(Sample)), ]




wq$Parkin <- wq$Sample
wq$Parkin <- gsub("[[:digit:]]", "", wq$Parkin)
wq$Parkin <- gsub("^P", "Parkinson", wq$Parkin)
wq$Parkin <- gsub("^H", "Healthy", wq$Parkin)
we <- wq[-which(wq$mmol_per_gram_dry_weight_per_hour == 0),]
new_variables <- c("Phenylacetic acid", "Indole","L_Tryptophan", "Fructose", "Myristic acid","Acetate" , "Butyrate", "Propionate",  "3-Methyl-2-Oxovaleric Acid", "N-Acetylneuraminic acid")
ex_reac <- c("EX_pac(e)", "EX_indole(e)", "EX_trp_L(e)", "EX_fru(e)", "EX_ttdca(e)", "EX_ac(e)", "EX_but(e)", "EX_ppa(e)","EX_3mop(e)", "EX_acnam(e)")


#taken from: https://stackoverflow.com/questions/19424709/r-gsub-pattern-vector-and-replacement-vector
for (i in seq_along(ex_reac)){
  we$sub <- gsub(ex_reac[i], new_variables[i], we$sub, fixed = T)
}


resultado <- "C:/Users/57322/Desktop/cod/phy/"
rxn.test.dti <- data.table()

###############################################################################################
listaazucare <- vector("list") 
listaazucares <- vector("list")

for(i in new_variables){
  #i ="N-Acetylneuraminic acid"
  #i ="EX_trp_L(e)"
  print(i)
  
  hc <- we[which(we$sub == i),]
  
  
  
  hc1 <- aggregate(mmol_per_gram_dry_weight_per_hour~spec+time+sub+Parkin,hc,mean)
  Healthys<-filter(hc1, Parkin == "Healthy")
  Parkinsons<-filter(hc1, Parkin == "Parkinson")
  
  hc2_h_cero <-Healthys[Healthys$mmol_per_gram_dry_weight_per_hour >= 0,]
  if (nrow(hc2_h_cero)!=0){ 
    hc2_h_tail <-head(hc2_h_cero,10) }else {
      hc2_h_tail <-hc2_h_cero
     
    } 
  
  
  
  hc2_h_ceromas <-Healthys[Healthys$mmol_per_gram_dry_weight_per_hour <= 0,]
  if (nrow(hc2_h_ceromas)!=0){ 
    hc2_h_tails <-head(hc2_h_ceromas,10)}else {
      hc2_h_tails <-hc2_h_ceromas
    } 
  
  
  hc2_p_cero <-Parkinsons[Parkinsons$mmol_per_gram_dry_weight_per_hour >= 0,]
  if (nrow(hc2_p_cero)!=0){ 
    hc2_p_tail <-head(hc2_p_cero,10) }else {
      hc2_p_tail <-hc2_p_cero
    } 
  
  
  
  hc2_p_ceromas <-Parkinsons[Parkinsons$mmol_per_gram_dry_weight_per_hour <= 0,]
  if (nrow(hc2_p_ceromas)!=0){ 
    hc2_p_tails <-head(hc2_p_ceromas,10)}else {
      hc2_p_tails <-hc2_p_ceromas
    } 
  
  
  
  ##unir
  
  
  
  
  #
  exists("hc2_h_tail") && is.data.frame(get("hc2_h_tail"))
  
  if (exists("hc2_h_tail")){
    hc2_t <- rbind(hc2_h_tail, hc2_p_tail)
    
    #1 tail
    zt<-hc2_t %>%
      group_by(Parkin) %>%
      mutate(percent = (abs(mmol_per_gram_dry_weight_per_hour)/sum(abs(mmol_per_gram_dry_weight_per_hour)))*100)
    
    
    #, size = Parkin
    afg <-ggplot(zt, aes(x = Parkin, y = percent, color = Parkin, fill = spec, size = Parkin)) +
      geom_col(size=1) +
      geom_text(aes(label =  paste(round(percent,2),"%")), size = 4, position = position_stack(vjust = 0.5), fontface = "bold",colour = "white", vjust=0.5) +
      coord_polar(theta = "y") +
      scale_y_continuous(breaks = NULL) +
      scale_color_viridis(discrete = TRUE, option = "G")+
      scale_colour_viridis_d(option = "plasma")+
      labs(x = NULL, y = NULL) +
      theme_void()
    
    afg <- afg+ labs(color = "Phenotype")#, fill = "Species")
    afg <- afg+ labs(fill = "Organism")
    afg <- afg+
      labs(title = i)
    afg <- afg+ theme(axis.text.x = element_text(face="plain", color="black", 
                                                 size=16, angle=0),legend.title=element_text(size=16), legend.text=element_text(face="italic",size=16),)
    #afg1 <- afg1+ theme_classic( base_size = 8 ) 
    
    ggsave(afg,                                        #nombre de la gráfica en R
           file=paste(resultado,"dona_posit",i, ".jpeg", sep=''), #Nombre del jpeg
           height=20, width=30, units="in", dpi=300, bg = "white", limitsize = FALSE)
    listaazucare[[i]] <- afg  
  } else {}
  #correr 2
  

  exists("hc2_h_tails") && is.data.frame(get("hc2_h_tails"))
  if (exists("hc2_h_tails")&& is.data.frame(get("hc2_h_tails"))){
    hc2_h <-  rbind(hc2_h_tails, hc2_p_tails)
    
    zt1<-hc2_h %>%
      group_by(Parkin) %>%
      mutate(percent = (abs(mmol_per_gram_dry_weight_per_hour)/sum(abs(mmol_per_gram_dry_weight_per_hour)))*100)
    #, size = Parkin....+
    
    
    afg1 <-ggplot(zt1, aes(x = Parkin, y = percent,   color = Parkin, fill = spec, size = Parkin)) +
      geom_col(size=1) +
      geom_text(aes(label =  paste(round(percent,2),"%")), size = 4, position = position_stack(vjust = 0.5), fontface = "bold",colour = "white", vjust=0.5) +
      coord_polar(theta = "y") +
      scale_y_continuous(breaks = NULL) +
      #annotate(geom = "text", x=0.5, y=0, label = Parkin, size = 16, color = "grey") +
      scale_color_viridis(discrete = TRUE, option = "G")+
      scale_colour_viridis_d(option = "plasma")+
      labs(x = NULL, y = NULL) +
      theme_void()
    
    ##############
    afg1 <- afg1+ labs(color = "Phenotype")#, fill = "Species")
    afg1 <- afg1+ labs(fill = "Organism")
    afg1 <- afg1+
      labs(title = i)
    afg1 <- afg1+ theme(axis.text.x = element_text(face="plain", color="black", 
                                                   size=16, angle=0),legend.title=element_text(size=16), legend.text=element_text(face="italic",size=16),)
    #afg1 <- afg1+ theme_classic( base_size = 8 ) 
    ggsave(afg1,                                        
           file=paste(resultado,"donas",i, ".jpeg", sep=''), 
           height=20, width=30, units="in", dpi=300, bg = "white", limitsize = FALSE)
    listaazucares[[i]] <- afg1  
  } else {}
  
  #listaazucares[[i]] <- afg1  
  datavector<-c("zt","zt1","hc2_t", "hc2_h","afg1","afg","hc2_h_tail", "hc2_p_tail", "hc2_h_tails", "hc2_p_tails")
  for(j in datavector){
    if (exists(j)){
      rm(j)
    }  else {}
    
    
    
  }
}
###fin
#este
#Positives-production
h1<- as.grob(listaazucare[[1]]+ guides(color="none")) 
h2<- as.grob(listaazucare[[2]]+ guides(color="none"))
h3<- as.grob(listaazucare[[3]]+ guides(color="none")) 
h4<- as.grob(listaazucare[[4]]+ guides(color="none"))
h5<- as.grob(listaazucare[[5]]+ guides(color="none"))
h6<- as.grob(listaazucare[[9]])
h7<- as.grob(listaazucare[[10]])


layout <- rbind(c(1, 2),
                c(4, 5),
                c(6,6))

milo<-grid.arrange(h1,h2, h4, h5,h6,layout_matrix=layout,  top = textGrob("Bacterial Metabolites production",gp=gpar(fontsize=20)))
ggsave("C:/Users/57322/Desktop/cod/donut_metabol_pos.png", plot = milo, width = 25, height = 15, units = "in", dpi = 300)


h8<- as.grob(listaazucare[[6]]+ guides(color="none"))
h9<- as.grob(listaazucare[[7]]+ guides(color="none"))
h10<- as.grob(listaazucare[[8]])


layout <- rbind(c(1, 2),
                c(3,3))

milo2<-grid.arrange(h8,h9,h10,layout_matrix=layout,  top = textGrob("Bacterial SCFAs production",gp=gpar(fontsize=20)))

ggsave("C:/Users/57322/Desktop/cod/donut_scfas_pos.png", plot = milo2, width = 25, height = 15, units = "in", dpi = 300)

#negatives-uptake
hs1<- as.grob(listaazucares[[1]]+ guides(color="none")) 
hs2<- as.grob(listaazucares[[2]]+ guides(color="none"))
hs3<- as.grob(listaazucares[[3]]+ guides(color="none")) 
hs4<- as.grob(listaazucares[[4]]+ guides(color="none"))
hs5<- as.grob(listaazucares[[5]]+ guides(color="none"))
hs6<- as.grob(listaazucares[[9]]+ guides(color="none"))
hs7<- as.grob(listaazucares[[10]])


layout <- rbind(c(2, 3),
                c(4, 5),
                c(6, 7))

moly<-grid.arrange(hs2,hs3,hs4, hs5,hs6, hs7,layout_matrix=layout,  top = textGrob("Bacterial Metabolites uptake",gp=gpar(fontsize=20)))

ggsave("C:/Users/57322/Desktop/cod/donut_metabol_neg.png", plot = moly, width = 25, height = 15, units = "in", dpi = 300)

hs8<- as.grob(listaazucares[[6]]+ guides(color="none"))
hs9<- as.grob(listaazucares[[7]]+ guides(color="none"))
hs10<- as.grob(listaazucares[[8]])


layout <- rbind(c(1, 2),
                c(3, 3))

moli2<-grid.arrange(hs8,hs9,hs10,layout_matrix=layout,  top = textGrob("Bacterial SCFAs uptake",gp=gpar(fontsize=20)))


ggsave("C:/Users/57322/Desktop/cod/donut_scfas_neg.png", plot = moli2, width = 25, height = 15, units = "in", dpi = 300)

##############################################################################################
