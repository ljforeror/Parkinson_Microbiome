
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

####################################### Data bases

data <- read.table("F:/articulo/decompri/data.txt", header=T, row.names=1, check.names=T)
metadata <- read.table("F:/articulo/decompri/Metadata_english.txt")
metadatas <- read.table("F:/articulo/decompri/metadatas.txt",header=T, row.names=1,check.names=F, sep=";")
data_sincorreccion_kilocal <- read.csv("F:/articulo/decompri/data_sincorreccion_kilocal.txt", sep="")
otu_mat<- read.table("F:/articulo/decompri/_joha_ZOTUstab_la_1.txt", header=T, row.names=1,check.names=F, sep="\t")
tax_mat <-read.table("F:/articulo/decompri/_joha_taxo_2022.csv",  header=T, row.names=1,check.names=F, sep=",")
samples_JF <- read.table("F:/articulo/decompri/samples_metadata_english.txt")
tree <- read.tree("F:/articulo/decompri/_joha_clusters.tree")



########################################## Demographic analysis

################Fisher test
# choosing the important variables
borrar <- c("samples","colors","ID", "Fiber","caffeine", "PC1","PC2","PC1_d","PC2_d", "number")
metasi <- metadata[ , !(names(metadata) %in% borrar)]
fishersitas <- as.character(c("Sex" ,"Smoke", "Alcohol", "levodopacarvidopa","Pramipexole", "Phenobarbital", "Rotigotine","Mirapex","Metformin", "Level_PD"))
fisher_variables <- data.table()
table_all <- vector(mode = "list",length = length(fishersitas))
idx.con <- grep("^c",rownames(metasi))
idx.pat <- grep("^p",rownames(metasi))

##organizing table

g.pat <- metadata[,"phenotype"]
phenotype <- metadata$phenotype
for(i in fishersitas){
  
  
  asi <- table(phenotype,as.factor(metasi[[i]]) )
  
  ese <- prop.table(asi)
  ise <- ese*100
  table_alles <- cbind(asi,ise)
  table_all[[i]] <- as.data.frame(t(table_alles[,4]))
  si <- as.data.frame(t(table_alles[,4]))
  colnames(si)<- c("Healthy%", "Patients%")
  num <- as.data.frame(t(asi[,2]))
  
  fisher <- fisher.test(asi)
  
  
  
  fisher_variables <- rbind(fisher_variables, data.table(Variables= i , num, si, fisher=fisher$p.value))
  
}
fisher_variables[,fisher.adjust:=p.adjust(fisher)]
fisher_variables
fisher_variables<- fisher_variables %>% mutate_if(is.numeric, ~round(., 3))

#saving
write.table(fisher_variables, file ="C:/Users/57322/Desktop/cod/demografic_fisher.csv", dec=",")

filter_table_correc_non <- fisher_variables[fisher<=0.05]
#correction
filter_table_correc_noni <- fisher_variables[fisher.adjust<=0.05]

sie<- table(phenotype,metasi[, "Level_PD"]) 
fisher <- fisher.test(sie)

#correction
pairwise_fisher_test(as.matrix(sie), p.adjust.method = "fdr")

######################wilcoxon Test

#organizing table and choosing variables
wc_variables <- data.table()
wilcoxonsitas <-c("Stoolsweek","Scale", "WEBSTER", "Weight", "Calf_Perimeter")

for(i in wilcoxonsitas){
  #print(i)
  
  g.con <- as.numeric(metasi[idx.con,i])
  g.pat <- as.numeric(metasi[idx.pat,i])
  
  wc <- wilcox.test(g.con,g.pat)
  ks <- ks.test(g.con, g.pat)
  g.con <- as.data.frame(as.numeric(metasi[idx.con,i]))
  g.pat <- as.data.frame(as.numeric(metasi[idx.pat,i]))
  median_PD <- get_summary_stats(g.pat, type="median_iqr")
  median_controls <- get_summary_stats(g.con, type="median_iqr")
  wc_variables <- rbind(wc_variables, data.table(Variables=i, p.value=wc$p.value,Median_PD = median_PD, Median_controls= median_controls ))
  
}
adf<-wc_variables[,wc.adjust:=p.adjust(p.value)]
adf<- adf %>% mutate_if(is.numeric, ~round(., 3))
#wc_variables
write.table(adf, file ="C:/Users/57322/Desktop/cod/demografic_wc.csv", dec=",")
#correction
filter_table_correc <- wc_variables[wc.adjust<=0.05]
filter_table_correc <- filter_table_correc[,-c(3,7)]
filter_table_correc
#without correction
filter_table_correc_non <- wc_variables[p.value<=0.05]
filter_table_correc_non <- filter_table_correc_non[,-c(3,7)]
filter_table_correc_non

###############END

################################################ Dietary analysis
## Descriptive Statistics
### Wilcoxon test
#### Significant differences between patients and healthy
##### Wilcoxon test noncorrected metabolites from the diet

#organizing table and choosing important variables
data$samples<- NULL
data$VITAMINE.D..MGR.<- NULL
patients <- data[,"Phenotype"]==0
healthy <- data[,"Phenotype"]==1

sigMat <- matrix(0,ncol(data)-1,3)
colnames(sigMat) <- c("Median healthy","Median patients","p.value")
rownames(sigMat) <- colnames(data)[-(1)]

for (i in 2:ncol(data)){
  
  valPatients <- data[patients,i]
  valHealthy <- data[healthy,i]
  sigMat[i-1,1] <- median(valHealthy)
  sigMat[i-1,2] <- median(valPatients)
  pval <- (wilcox.test(valPatients,valHealthy))$p.value
  sigMat[i-1,3] <- pval
}
a <- as.data.frame(sigMat)
tabli <- a[which(a$p.value < 0.05), ]
tabli

kable(tabli, escape = F, caption = "Wilcoxon_test_non_corrected_metabolites_diet")%>%
  kable_styling(latex_options = "hold_postition")

write.table(format(tabli,scientific=FALSE, digits=4), file ="C:/Users/57322/Desktop/cod/Wilcoxon_test_non_corrected_metabolites_diet.csv", dec=",")

##### Wilcoxon test FDR corrected metabolites from the diet

sigMat <- matrix(0,ncol(data)-1,3)
colnames(sigMat) <- c("Median healthy","Median patients","padj")
rownames(sigMat) <- colnames(data)[-(1)]

for (i in 2:ncol(data)){
  
  valPatients <- data[patients,i]
  valHealthy <- data[healthy,i]
  sigMat[i-1,1] <- median(valHealthy)
  sigMat[i-1,2] <- median(valPatients)
  
  pval <- (wilcox.test(valPatients,valHealthy))$p.value
  
  sigMat[i-1,3] <- pval
  
}

sigMat[,3] <- p.adjust(sigMat[,3],method="fdr")
a <- as.data.frame(sigMat)
tabli_2 <- a[which(a$padj < 0.05), ]
tabli_2
write.table(format(tabli_2,scientific=FALSE, digits=4), file ="C:/Users/57322/Desktop/cod/Wilcoxon_test_corrected_metabolites_diet.csv", dec=",")


##################### Graphics

data$Phenotype <- gsub("1", "Healthy", data$Phenotype)
data$Phenotype <- gsub("0", "Patient", data$Phenotype)

###Long format

datas <- data
datas$ID <- rownames(data)


diet_met <- gather(datas,
                   key = "Metabolites_diet",
                   value = "Uptake_Amount",
                   -c("ID", "Phenotype"))

new_diet_met <- diet_met %>% 
  select(Phenotype, ID, Uptake_Amount, Metabolites_diet ) %>%
  filter(Metabolites_diet %in% c("Trans_fat", "Carbohydrates", "Potassium"))
new_diet_met$Phenotype <- gsub("1", "Healthy", new_diet_met$Phenotype)
new_diet_met$Phenotype <- gsub("0", "Patient", new_diet_met$Phenotype)

p <- ggplot(data = new_diet_met, aes(x=Phenotype, y=Uptake_Amount, fill=Phenotype)) +   geom_boxplot() + 
  geom_jitter(width=0.1,alpha=0.2) +
  xlab("Phenotype")+ 
  facet_wrap(vars(Metabolites_diet), scales = "free_y") +
  scale_fill_manual(values=c("#CD6889", "#8B668B"))+
  theme_classic()

p
#ggsave("C:/Users/57322/Desktop/cod/Uptake_Amount_300.png", plot = p, width = 8, height = 5, units = "in", dpi = 300)
#dev.off()


###################### PCA

Phenotype <- as.vector(data$Phenotype)
datas$Phenotype <- gsub("1", "Healthy", datas$Phenotype)
datas$Phenotype <- gsub("0", "Patient", datas$Phenotype)

res.pca_1 <- prcomp(datas[,2:60], scale = TRUE)
fviz_eig(res.pca_1)
x <-fviz_pca_ind(res.pca_1, axes = c(1, 2), label="none",
                 col.ind = datas$Phenotype, # Color by the quality of representation
                 
                 repel = TRUE,
                 addEllipses=TRUE, ellipse.level=0.95, palette = c("#CD6889", "#8B668B")
)   # Avoid text overlapping
pca_dieta <- fviz_pca_ind(res.pca_1,
                          col.ind = data$Phenotype, # Color by the quality of representation
                          #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                          repel = TRUE, addEllipses=TRUE  )   # Avoid text overlapping

y<-fviz_contrib(res.pca_1, choice = "var", axes = 1, top = 15)
z <- fviz_contrib(res.pca_1, choice = "var", axes = 2, top = 15)
zz <- fviz_pca_biplot(res.pca_1, axes = c(1, 2), addEllipses = TRUE, repel = TRUE, select.var = list(contrib = 20),  label="var",
                      col.var = "#63B8FF", # Variables color
                      col.ind = datas$Phenotype,  # Individuals color
                      palette = c("#CD6889", "#8B668B")
)
fviz_pca_biplot(res.pca_1,axes = c(2, 3), repel = TRUE, select.var = list(contrib = 20),  label="var",
                col.var = "#CD6889", # Variables color
                col.ind = "#8B668B"  # Individuals color
)

#####Save information
#ggsave("C:/Users/57322/Desktop/cod/pca_nutrients_1_300.png", plot = x, width = 8, height = 5, units = "in", dpi = 300, bg = "white")
#dev.off()
#ggsave("C:/Users/57322/Desktop/cod/pca_nutrients_2_300.png", plot = y, width = 8, height = 5, units = "in", dpi = 300, bg = "white")
#dev.off()
#ggsave("C:/Users/57322/Desktop/cod/pca_nutrients_3_300.png", plot = z, width = 8, height = 5, units = "in", dpi = 300, bg = "white")
#ggsave("C:/Users/57322/Desktop/cod/pca_nutrients_4_300.png", plot = zz, width = 8, height = 5, units = "in", dpi = 300, bg = "white")
#colors <- c("#CD6889", "#8B668B")


#arrange <- ggarrange(zz, new_diet_met,ord3, ord4,  ncol = 2, nrow = 2, labels  = c("A","B", "C", "D"),  widths = c(0.5, 0.50))
#ggsave("C:/Users/57322/Desktop/cod/arrange_diet_300.png", plot = arrange, width = 20, height = 12, units = "in", dpi = 300)



#Open files: data = metabolite information, metadata = clinical and anthropometric variables
# Convert titles from upper case to lower case (all tables)

#colnames(data)=tolower(colnames(data))
#colnames(metadata)=tolower(colnames(metadata))

#data[,1]=tolower(data[,1])
#data[1]
data$samples <- rownames(data)
#colnames(data)[1] <- "Samples"

metadata$phenotype <- NULL
metadata$colors <- NULL
metadata$id <- NULL

#Merge of metadata and diet info

datosMERGE=merge(metadata,data,by="samples")

#Transfrom to numeric variables

datosMERGE$Sex=ifelse(datosMERGE$Sex=="M",0,1)
datosMERGE$Smoke=ifelse(datosMERGE$Smoke=="NO",0,1)
datosMERGE$Alcohol=ifelse(datosMERGE$Alcohol=="NO",0,1)

datosMERGE$Level_PD[datosMERGE$Level_PD=="No"]=0
datosMERGE$Level_PD[datosMERGE$Level_PD=="Early"]=1
datosMERGE$Level_PD[datosMERGE$Level_PD=="Moderate"]=2
datosMERGE$Level_PD[datosMERGE$Level_PD=="Severe"]=3

as.factor(datosMERGE$Level_PD)

#organizing important varibles
rownames(datosMERGE)=datosMERGE$samples
datosMERGE=datosMERGE[,-1]

datosMERGE$Sex=as.numeric(datosMERGE$Sex)
datosMERGE$Level_PD=as.numeric(datosMERGE$Level_PD)

datosMERGE=datosMERGE[,-c(1,9,11:39)]
datosMERGE$number <- NULL
datosMERGE$todo <- NULL
datosMERGE$micro_16sinson <- NULL
apply(datosMERGE,2,sd)

summary(datosMERGE)


metadatas <- datosMERGE
rownames(metadatas) <- c("h1","h10","h11","h12","h13","h14","h15","h16","h17","h18","h19","h2","h20","h21","h22","h23","h24","h25","h3","h4","h5","h6","h7","h8","h9","p1","p10","p11","p12","p13","p14","p15","p16","p17","p18","p19","p2","p20","p21","p22","p23","p24","p25","p3","p4","p5","p6","p7","p8","p9")
metadatas$PC1 <- NULL
metadatas$PC2 <- NULL
metadatas$PC1_d <- NULL
metadatas$PC2_d <- NULL
metadatas$kilocal <- NULL
metadatas$Fiber <- NULL
metadatas$caffeine <- NULL
metadatas$Phenotype <- NULL
metadatas$COPD <- NULL
metadatas$Vitamin_A.1 <- NULL
metadatas$Sex <- NULL
metadatas$Smoke  <- NULL
metadatas$Alcohol <- NULL



#Correlation and p.value

class(metadatas)
setDT(metadatas)

correlacion <- tab_corr(metadatas)
tab_corr(metadatas,p.numeric = TRUE, corr.method = "spearman")

tab_corr(metadatas,
         corr.method = "spearman",
         p.numeric = TRUE,
         file = "Correlacion_Matrix.doc")
getwd()


#adjust table by p.values

p_ajustado <- rcorr.adjust(as.matrix(metadatas), type = "spearman")
valores_p <- p_ajustado$P
#round(valores_p,3)
#View(valores_p)
class(valores_p)


corr_final <- p_ajustado$R
correlaciones <-corr_final$r
filtro_corr <- correlaciones[1:8,9:66]
final_corr <- round(filtro_corr,3)

row.names(final_corr) <- c("Age", "Weight", "Stools/week", "Bristol_score","Disease_years","Webster_Score","Level_PD","Calf_Perimeter")

#View(final_corr)
#####Save information
#ggsave("figura_johanna.jpg", plot = (imagen), width = 20, height = 12, units = "in", dpi = 300)
#ggsave("imagen2_johanna.png", dpi = 300)



plot(valores_p)
filtro <- valores_p[1:8,9:66]


#################
#####Save information
#Filtrar tabla de correlacion para que solo tenga valores significativos p
#se seleccionan los valores p significativos a partir de la siguiente tabla
write.table(filtro, "filtro.csv", sep = ';', dec = '.')

# then we build the correlation table and eliminate the columns so that only the ones we chose based on the p-value remain
write.table(format(correlaciones, scientific=FALSE,digits=4), "tabla_correlaciones.csv", sep = ';', dec = '.')
write.table(format(final_corr, scientific=FALSE,digits=4), "final_corr.csv", sep = ';', dec = '.')

#Nutrients vs Nutrients Only
nutrientes_p <- valores_p[9:66,9:66]
#write.table(format(nutrientes_p, scientific=FALSE,digits=4), "nutrientesvsnutrientes.csv", sep = ';', dec = '.')

corclinicasdieta <- cor(metadatas, method = "spearman")

#heatmap
heatmap(final_corr, scale="column", col = heat.colors(256))


##PCA Biplot diet Healthy vs Patient 
col<- colorRampPalette(c("red", "white", "blue"))(10)
transform_p <- function(x) {
  y <- 0.91 - (0.82) * (1 - exp(-3.82 * x))
  y
}

#corplot
imagen <- corrplot(final_corr,
                   method = "circle",
                   type = "full",
                   tl.pos = "tl",
                   tl.col = "black",
                   order = "original",
                   tl.cex = 0.6,
                   col = col,
                   size_vector = transform_p(final_corr),
                   cl.cex= 0.8)

### Variables at 45grades
p1 <- corrplot(final_corr,
               method = "circle",
               type = "full",
               tl.pos = "tl",
               tl.col = "black",
               tl.srt = 45,
               order = "original",
               tl.cex = 0.9,
               col = col,
               cl.cex= 0.9,
               size_vector = transform_p(final_corr))



grid.echo()
P1 <- grid.grab()

grid.draw(P1) 
P1 <- editGrob(P1,
               gPath("background"), grep = TRUE,
               gp = gpar(fill = NA))

#####Save information
#ggsave(filename = "C:/Users/57322/Desktop/cod/p1.png", plot = replayPlot(p1), dpi = 300)
#ggsave("C:/Users/57322/Desktop/cod/Correlation_5matrix.png",
#replayPlot(p1),
#width = 3.5,
#height = 3,
#dpi = 300
#)



#####
#file_path= "C:/Users/57322/Desktop/cod/Correlation_matrix.png"
#png(height=1000, width=3500, file=file_path, type = "cairo", res = 300, units = "px", pointsize = 12)
##different format corrplot

p1 <- {corrplot(final_corr,
                method = "circle",
                type = "full",
                tl.pos = "tl",
                tl.col = "black",
                tl.srt = 45,
                order = "original",
                tl.cex = 0.8,
                col = col,
                cl.cex= 0.8,
                size_vector = transform_p(final_corr))
  recordPlot()}
# Then
dev.off()


#Loading required package: grid
myPlot1 = p+ theme(legend.position="none", strip.text = element_text(size = 16, colour = "black")) + font("legend.text", size = 16, color = "black", face = "bold")+ font("axis.title", size = 18, color = "black", face = "bold") + font("axis.text", size = 16, color = "black", face = "bold")
myPlot2=zz+ theme(strip.text = element_text(size = 16, colour = "black")) + font("legend.text", size = 16, color = "black", face = "bold")+ font("axis.title", size = 18, color = "black", face = "bold") + font("axis.text", size = 16, color = "black", face = "bold")

myPlot=list(P1)

#### Arrange final of the section
re<-ggdraw() +
  draw_plot(P1, x = 0, y = 0, width = 1, height = .7) +
  draw_plot(myPlot2, x = .5, y = .5, width = .5, height = .4) +
  draw_plot(myPlot1, x = 0, y = 0.5, width = 0.5, height = 0.4) +
  draw_plot_label(label = c("A", "C", "B"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))

#####Save information
#ggsave("C:/Users/57322/Desktop/cod/arrange_2_final_300.png", plot = re, width = 15, height = 12, units = "in", dpi = 300)

#############################################################################################
############################################### 16S Taxonomic analysis
## some parts are taken from : https://joey711.github.io/phyloseq/

##organizing samples
df_ordenado_otu <- otu_mat[,order(colnames(otu_mat)) ]
#merge metadatas
samples_JF <- merge(samples_JF, metadata)

#check if the names of my tables are correct
rownames(samples_JF) <- samples_JF$samples
colnames(otu_mat) <- samples_JF$samples
all(rownames(samples_JF) == colnames(otu_mat))
colnames(otu_mat)
colnames(samples_JF)
colnames(tax_mat)

#convet to matrix
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#put Phenotype variable a factor
colnames(samples_JF)[1] <- "Phenotype"
Phenotype <- samples_JF$Phenotype

# phyloseq-object
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_JF)

#mostrar colores
samples_JF$colors

# phyloseq
objphil <- phyloseq(OTU, TAX, samples, tree)

#without rarefy
plot_heatmap(objphil, taxa.label="Phylum")

#  color column to be of type "character" for plots
samples_JF$colors<- as.character(samples_JF$colors)

#Remove low occurence / abundance zOTU more than 10 sequences in total and appearing in more than 1 sample

objphil <- filter_taxa(objphil, function(x) sum(x >= 10) > (2), prune =  TRUE)

#RAREFACCION
vegan::rarecurve(t(otu_table(objphil)), step=100, col = sample_data(objphil)$colors, lwd=2, ylab="zOTUs", label=F)
abline(v=(min(rowSums(t(otu_table(objphil))))))
###prevalence more reads
#Remove taxa not seen more than 3 times in at least 20% of the samples.
balea_2 <-filter_taxa(objphil, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

###get table from balea_2
tabla_zotus_filtrada <- balea_2@otu_table
colSums(tabla_zotus_filtrada)
mean(colSums(tabla_zotus_filtrada))
mean(colSums(otu_mat))
#saving info
#write.table(format(tabla_zotus_filtrada,scientific=FALSE, digits=4), file ="C:/Users/57322/Desktop/cod/filtrada_balea_2.csv", dec=",")



#####################################normalization
#taken from: https://micca.readthedocs.io/en/1.7.1/phyloseq.html
# 90% of the abundance of the sample with less reads
ps.rarefied2 <- phyloseq::rarefy_even_depth(balea_2, rngseed=123, sample.size=0.9*min(sample_sums(balea_2)), replace=F)
ps.rarefied <- ps.rarefied2

vegan::rarecurve(t(otu_table(ps.rarefied)), step=200, col = sample_data(ps.rarefied)$colors, lwd=2, ylab="ASVs", label=F)
abline(v=(min(rowSums(t(otu_table(ps.rarefied))))))

vegan::rarecurve(t(otu_table(ps.rarefied)[1:400,]), step = 20,  col =
                   sample_data(ps.rarefied)$colors, ylab="zOTUs", cex =      0.6,
                 main = "Rarecurve")
rare <- vegan::rarecurve(t(otu_table(ps.rarefied)[1:400,]), step = 20,  col =
                           sample_data(ps.rarefied)$colors, ylab="zOTUs", cex =      0.6,
                         main = "Rarecurve")
#saving info
#ggsave("C:/Users/57322/Desktop/cod/rare_2_300.png", plot = rare, width = 20, height = 12, units = "in", dpi = 300)

#Relative abundance
mean(colSums(ps.rarefied@otu_table))
#Get mean
ps_rel_abund = phyloseq::transform_sample_counts(balea_2, function(x){x / sum(x)})
#####taxa con 0.1 %

physeqrF = filter_taxa(ps_rel_abund, function(x) mean(x) > 0.01,TRUE)

#heatmap mas de 0.1
carbom_abund <- phyloseq::filter_taxa(ps_rel_abund, function(x) sum(x > 0.10) > 0, TRUE)
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))((25))
plot_heatmap(carbom_abund, method = "NMDS", distance = "bray") + facet_wrap(~Phenotype, scales = "free_x", nrow = 1)




fantaxtic_bar(carbom_abund, color_by = "Family", label_by = "Genus", facet_by = NULL, grid_by = NULL, other_color = "Grey") -> ptop15

fantaxtic_bar(carbom_abund, color_by = "Family", label_by = "Genus", facet_by = NULL, grid_by = NULL, other_color = "Grey") -> ptop20


###abundances rarified
#taked from: https://github.com/joey711/phyloseq/issues/694
#average to normalized
ps_rel_abund1 = phyloseq::transform_sample_counts(ps.rarefied, function(x){x / sum(x)})



#heatmap more than 0.1 normalized
carbom_abund1 <- filter_taxa(ps_rel_abund1, function(x) sum(x > 0.05) > 0, TRUE)
carbom_abund2 <- filter_taxa(ps_rel_abund1, function(x) sum(x > 0.1) > 0, TRUE)

zotus_good_abu <-rownames(carbom_abund1@otu_table)

###save
sw<-stargazer(carbom_abund1@otu_table,                 # Export txt
              summary = FALSE,
              type = "text",
              out = "C:/Users/57322/Desktop/cod/carbom_abund1.txt")


plot_bar(ps_rel_abund1, fill="Phylum") + facet_wrap(~Phenotype, scales = "free_x", nrow = 1)+ theme_bw()
plot_bar(carbom_abund1, fill="Family") + facet_wrap(~Phenotype, scales = "free_x", nrow = 1)


ppk = plot_bar(carbom_abund1, "Family", fill="Genus", facet_grid=~Phenotype)
ppk + geom_point(aes(x=Family, y=Abundance), color="black", position="jitter", size=3)

plot_bar(carbom_abund1, "Family", fill="Genus", facet_grid=~Phenotype)
ppk
#########plot abundance relative
familia<-phyloseq::plot_bar(ps_rel_abund1, fill = "Family") +
  geom_bar(aes(color = Family, fill = Family), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Phenotype, scales = "free") +
  theme(panel.background = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), strip.text = element_text(size = 16, margin = margin()))
familia<-familia+ font("legend.title", size = 20, color = "black", face = "bold")+  font("legend.text", size = 14, color = "black", face = "italic")+ font("axis.title", size = 18, color = "black", face = "bold")+  font("axis.text", size = 15, color = "black", face = "italic") + font("subtitle", color = "black", face = "italic", size = 14)

Phylum<- phyloseq::plot_bar(ps_rel_abund1, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Phenotype, scales = "free") +
  theme(panel.background = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), strip.text = element_text(size = 16, margin = margin()))
Phylum<-Phylum+ font("legend.title", size = 20, color = "black", face = "bold")+  font("legend.text", size = 14, color = "black", face = "italic")+ font("axis.title", size = 18, color = "black", face = "bold")+  font("axis.text", size = 15, color = "black", face = "italic") + font("subtitle", color = "black", face = "italic", size = 14)


genus<- phyloseq::plot_bar(ps_rel_abund1, fill = "Genus") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Phenotype, scales = "free") +
  theme(panel.background = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#####saving plots
ggsave("C:/Users/57322/Desktop/cod/familia_300.png", plot = familia, width = 20, height = 12, units = "in", dpi = 300)
ggsave("C:/Users/57322/Desktop/cod/Phylum_300.png", plot = Phylum, width = 20, height = 12, units = "in", dpi = 300)
ggsave("C:/Users/57322/Desktop/cod/genus_300.png", plot = genus, width = 20, height = 12, units = "in", dpi = 300)
###################################################
#######other plot relative abundance
#based on https://github.com/joey711/phyloseq
ps1 <- merge_samples(objphil, "Phenotype")
ps1 <- transform_sample_counts(ps1, function (x) x / sum(x))

#Then aggregate zOTUs to the family level. Just to plot o illustration of zotus, or other taxa. Some info in metada is lost

ps2 <- tax_glom(ps1, "Family")


#keep families that have a proportion > 0.1 in at least one group

ps4 <- filter_taxa(ps2, function (x) max(x) > 0.1)

#wanted the top 15 families as defined by the mean proportion of the family across Group's.First list of the families you want to keep and use prune_taxa().

top_families <- names(sort(taxa_sums(ps2), decreasing = TRUE)[1:15])
ps5 <- prune_taxa(top_families, ps2)

#top_families will actually be a list of original zOTU names 
#check documentation to better implementation


df <- psmelt(ps5)
micro_16sinson <- as.factor(df$micro_16sinson)
colourCount =length(unique(df$Family))


qz<-ggplot(df, aes(Sample, Abundance, fill = Family)) +
  geom_bar(stat="identity")+
  ylab("Relative Abundance")+
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(colourCount))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
qz<- qz+ font("legend.title", size = 20, color = "black", face = "bold")+  font("legend.text", size = 14, color = "black", face = "bold")+ font("axis.title", size = 18, color = "black", face = "bold")+  font("axis.text", size = 15, color = "black", face = "bold") + font("subtitle", color = "black", face = "italic", size = 14)
qz
#save info
ggsave("C:/Users/57322/Desktop/cod/Relative_Abundance_Family_300.png", plot = qz, width = 20, height = 12, units = "in", dpi = 300)
#dev.off()
#######
########By Phylum
#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(ps_rel_abund1, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]

#Melt and plot
phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = Phenotype, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Phylum, scales = "free")+
  theme_classic()

#####By Family
#Agglomerate to phylum-level and rename
ps_Family <- phyloseq::tax_glom(ps_rel_abund, "Family")

#Melt and plot
phyloseq::psmelt(ps_Family) %>%
  ggplot(data = ., aes(x = Phenotype, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Family), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Family, scales = "free")+
  theme_classic()
##

psfa = subset_taxa(balea_2, Family == c("Streptococcaceae","Peptostreptococcaceae", "Verrucomicrobiaceae","Lactobacillaceae", "Lachnospiraceae"))
s <- phyloseq::psmelt(psfa) %>%
  ggplot(data = ., aes(x = Phenotype, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Family), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Family, scales = "free")+ scale_y_log10()+
  theme_classic()
s <- s + font("legend.title", size = 20, color = "black", face = "bold")+  font("legend.text", size = 14, color = "black", face = "italic")+ font("axis.title", size = 18, color = "black", face = "bold")+  font("axis.text", size = 15, color = "black", face = "italic") + font("subtitle", color = "black", face = "italic", size = 14)


ggsave("C:/Users/57322/Desktop/cod/Families_300.png", plot = s, width = 8, height = 8, units = "in", dpi = 300)
#dev.off()

psfa = subset_taxa(ps_rel_abund, Family == "Verrucomicrobiaceae")
verr<-phyloseq::psmelt(psfa) %>%
  ggplot(data = ., aes(x = Phenotype, y = Abundance, fill=Phenotype)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Genus), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Family, scales = "free") + scale_y_log10() +
  #scale_fill_manual(values=c("#CD6889", "#8B668B"))+
  theme_classic()
#arrange
arranget <- ggarrange(s + 
                        theme(legend.position="none", strip.text = element_text(size = 16, colour = "black")) + font("legend.text", size = 16, color = "black", face = "bold")+ font("axis.title", size = 18, color = "black", face = "bold") + font("axis.text", size = 16, color = "black", face = "bold"), 
                      
                      verr + 
                        font("legend.text", size = 14, color = "black", face = "bold")+ font("axis.title", size = 18, color = "black", face = "bold")+ font("axis.text", size = 15, color = "black", face = "bold")+
                        theme(legend.position="none")) 


draw_plot_label(label = c("A", "B", "C","D" ), size = 25,
                x = c(0, 0, 0.42,0.42), y = c( 1, 0.48, 1, 0.48)) # Add labels

#saving plot
ggsave("C:/Users/57322/Desktop/cod/box_fam_300.png", plot = arranget, width = 20, height = 12, units = "in", dpi = 300)
# Add labels
#################

#####################################Alpha Diversity
adiv2 <- data.frame(
  "Observed" = phyloseq::estimate_richness(ps.rarefied, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps.rarefied, measures = "Shannon"),
  "simpson" = phyloseq::estimate_richness(ps.rarefied, measures = "simpson"),
  "Phenotype" = phyloseq::sample_data(ps.rarefied)$Phenotype)
head(adiv2)


#Plot adiv measures
adiv2 %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson")) %>%
  ggplot(aes(x = Phenotype, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Phenotype), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none")

ay <- adiv2 %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson")) %>%
  ggplot(aes(x = Phenotype, y = value, color = Phenotype)) +
  geom_boxplot(outlier.color = NA) +
  scale_color_manual(values=c("#CD6889", "#8B668B")) +
  #scale_color_brewer(palette="Dark2") +
  geom_jitter(aes(color = Phenotype), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20))+ 
  theme_classic()
#saving info
ggsave("C:/Users/57322/Desktop/cod/alpha_300.png", plot = ay, width = 8, height = 8, units = "in", dpi = 300)
#dev.off()
####################################Significance Alpha diversity
#Wilcoxon test of location
wilcox.test(Observed ~ Phenotype, data = adiv2, exact = FALSE, conf.int = TRUE)

wilcox.test(Shannon ~ Phenotype, data = adiv2, conf.int = TRUE)   

wilcox.test(Simpson ~ Phenotype, data = adiv2, conf.int = TRUE)
## Beta Diversity

###############################Distances#######

# PCoA plot using the bray curtis as distance
braycur_dist <- phyloseq::distance(ps.rarefied, method="bray")
ordination = ordinate(ps.rarefied, method="PCoA", distance=braycur_dist)
plot_ordination(ps.rarefied, ordination, color="Phenotype") + theme(aspect.ratio=1)
ord <- plot_ordination(ps.rarefied, ordination, color="Phenotype") + theme(aspect.ratio=1)
ord + stat_ellipse(type = "norm", linetype =2)+ stat_ellipse(type = "t")+ theme_bw()


brq <-adonis2(braycur_dist~sample_data(ps.rarefied)$Phenotype, permutations = 9999)
br <-anosim(braycur_dist,sample_data(ps.rarefied)$Phenotype, permutations = 9999)

#option NMDS
ps_ordnmds <- ordinate(ps.rarefied, method = "NMDS", distance = braycur_dist)
plot_ordination(ps.rarefied, ps_ordnmds, color = "Phenotype")
ord2 <- plot_ordination(ps.rarefied, ps_ordnmds, color = "Phenotype")+ theme(aspect.ratio=1)
ord2 + stat_ellipse(type = "norm", linetype =2)+ stat_ellipse(type = "t")+ theme_bw()
#other alternative
ord2 <- plot_ordination(ps.rarefied, ps_ordnmds, color = "Phenotype", title="NMDS-Bray_curtis") + theme_bw()
ord2 + stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=Phenotype))
###

ord2 <- plot_ordination(ps.rarefied, ps_ordnmds, color = "Phenotype", title="NMDS-Bray_curtis") + 
  
  stat_ellipse(type="norm", alpha=0.4, linetype =2,  aes(fill="Phenotype"))  +
  stat_ellipse(type = "t")  +
  scale_fill_manual(breaks = c("Healthy", "Patients"),
                    values=c("#CD6889", "#8B668B")) +
  theme_classic()
gg_color_hue <- function(n){
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
color.names <- levels(ord2$data$Phenotype)
p4cols <- gg_color_hue(length(color.names))
names(p4cols) <- color.names
p4cols[1] <- "#CC0066"
p4cols[2] <- "#56B4E9"
q <- ord2 + scale_color_manual(values=p4cols)
#taken from: https://joey711.github.io/phyloseq/plot_ordination-examples.html

######################

#unifrac 
# NMDS plot using the weighted UniFrac as distance
wunifrac_dist = phyloseq::distance(ps.rarefied, method="unifrac", weighted=T)
ordination_wu = ordinate(ps.rarefied, method="NMDS", distance=wunifrac_dist)
plot_ordination(ps.rarefied, ordination_wu, color="Phenotype") + theme(aspect.ratio=1)
ord3 <- plot_ordination(ps.rarefied, ordination_wu, color="Phenotype", title="NMDS-Weight_Unifrac") + theme(aspect.ratio=1)
ord3 + stat_ellipse(type = "norm", linetype =2)+ stat_ellipse(type = "t")+ theme_bw()
############# improve plot
ord3 <- plot_ordination(ps.rarefied, ordination_wu, color = "Phenotype", title="NMDS-Weight_Unifrac")+ geom_point(size=3) +
  scale_color_manual(values = c("#CD6889", "#8B668B"))+
  stat_ellipse(type="norm", alpha=0.4, linetype =2,  aes(fill="Phenotype"))  +
  stat_ellipse(type = "t")  +
  scale_fill_manual(breaks = c("Healthy", "Patients"),
                    values=c("#CD6889", "#8B668B")) +
  theme_classic()+geom_vline(xintercept = 0, alpha=.25)+
  geom_hline(yintercept = 0, alpha=.25)
ord3 <- ord3 + theme(strip.text = element_text(size = 12, colour = "black", face = "italic"))+ font("axis.title", size = 18, color = "black", face = "bold") + font("legend.text", color = "black", face = "italic", size = 14) + font("axis.text", size = 15, color = "black", face = "bold") + font("legend.title", size = 20, color = "black", face = "bold") 


color.names <- levels(ord2$data$Phenotype)
p4cols <- gg_color_hue(length(color.names))
names(p4cols) <- color.names
p4cols[1] <- "#CD6889"
p4cols[2] <- "#8B668B"
e <- ord3 + scale_color_manual(values=p4cols)


##############
# NMDS plot using the unweighted UniFrac as distance
uwunifrac_dist = phyloseq::distance(ps.rarefied, method="unifrac", weighted=F)
ordination_uwu = ordinate(ps.rarefied, method="NMDS", distance=uwunifrac_dist)
plot_ordination(ps.rarefied, ordination_uwu, color="Phenotype") + theme(aspect.ratio=1)
ord4 <- plot_ordination(ps.rarefied, ordination_uwu, color="Phenotype", title="NMDS-UnWeight_Unifrac") + theme(aspect.ratio=1)
ord4 + stat_ellipse(type = "norm", linetype =2)+ stat_ellipse(type = "t")+ theme_bw()

############# improve plot
ord4 <- plot_ordination(ps.rarefied, ordination_uwu, color = "Phenotype", title="NMDS-UnWeight_Unifrac")+ geom_point(size=3) +
  scale_color_manual(values = c("#CD6889", "#8B668B"))+
  stat_ellipse(type="norm", alpha=0.4, linetype =2,  aes(fill="Phenotype"))  +
  stat_ellipse(type = "t")  +
  scale_fill_manual(breaks = c("Healthy", "Patients"),
                    values=c("#CD6889", "#8B668B")) +
  theme_classic()+geom_vline(xintercept = 0, alpha=.25)+
  geom_hline(yintercept = 0, alpha=.25)

ord4 <- ord4 + theme(strip.text = element_text(size = 12, colour = "black", face = "italic"))+ font("axis.title", size = 18, color = "black", face = "bold") + font("legend.text", color = "black", face = "italic", size = 14) + font("axis.text", size = 15, color = "black", face = "bold") + font("legend.title", size = 20, color = "black", face = "bold")

color.names <- levels(ord2$data$Phenotype)
p4cols <- gg_color_hue(length(color.names))
names(p4cols) <- color.names
p4cols[1] <- "#CD6889"
p4cols[2] <- "#8B668B"
l <- ord4 + scale_color_manual(values=p4cols)
#
#######Adonis other variables individually

ak1<-adonis2(wunifrac_dist~sample_data(ps.rarefied)$Phenotype +sample_data(ps.rarefied)$Age, permutations = 9999)
ak2<-adonis2(uwunifrac_dist~sample_data(ps.rarefied)$Phenotype +sample_data(ps.rarefied)$levodopacarvidopa, permutations = 9999)
ak3<-adonis2(uwunifrac_dist~sample_data(ps.rarefied)$Phenotype +sample_data(ps.rarefied)$WEBSTER, permutations = 9999)
ak4<-adonis2(uwunifrac_dist~sample_data(ps.rarefied)$Phenotype +sample_data(ps.rarefied)$levodopacarvidopa, permutations = 9999)
ak5<-adonis2(uwunifrac_dist~sample_data(ps.rarefied)$Phenotype +sample_data(ps.rarefied)$PC1_d, permutations = 9999)
ak6<-adonis2(uwunifrac_dist~sample_data(ps.rarefied)$Phenotype +sample_data(ps.rarefied)$disease_years, permutations = 9999)
ak7<-adonis2(uwunifrac_dist~sample_data(ps.rarefied)$Phenotype +sample_data(ps.rarefied)$Age, permutations = 9999)
ak8<-adonis2(uwunifrac_dist~sample_data(ps.rarefied)$Phenotype +sample_data(ps.rarefied)$Calf_Perimeter, permutations = 9999)
ak<-rbind(ak1[1:2,3:5],ak2[1:2,3:5],ak3[1:2,3:5],ak4[1:2,3:5],ak5[1:2,3:5],ak6[1:2,3:5],ak7[1:2,3:5], ak8[1:2,3:5])#, ak9[1:2,3:5])
write.table(format(ak,scientific=FALSE, digits=4), file ="C:/Users/57322/Desktop/cod/table_adonis.csv", dec=",")

############Anosim
uni <-anosim(uwunifrac_dist,sample_data(ps.rarefied)$Phenotype, permutations = 9999)
aswe<-anosim(wunifrac_dist,sample_data(ps.rarefied)$Phenotype, permutations = 9999)


#Adonis
# make a data frame from the sample_data
rrsampledf1 <- data.frame(sample_data(ps.rarefied))

# Adonis test
adonis(braycur_dist ~ Phenotype, data = rrsampledf1, permutations = 1000)
variable_group = get_variable(rrsampledf1, "Phenotype")
variable_group = anosim(uwunifrac_dist, variable_group)

uww <- adonis(uwunifrac_dist ~ sample_data(ps.rarefied)$Phenotype, data = rrsampledf1, permutations = 9999)
variable_group = get_variable(rrsampledf1, "Phenotype")
variable_group = anosim(uwunifrac_dist, variable_group)

wun<-adonis(wunifrac_dist ~ sample_data(ps.rarefied)$Phenotype, data = rrsampledf1, permutations = 9999)
variable_group = get_variable(rrsampledf1, "Phenotype")
variable_group = anosim(wunifrac_dist, variable_group)



plot((variable_group)
     ,main="ANOSIM wunifrac "
     ,las=1)
#########Other test
#anova
#beta
#Permutest


#wunifrac_dist
anova(betadisper(wunifrac_dist, rrsampledf1$Phenotype))
beta <- betadisper(wunifrac_dist, rrsampledf1$Phenotype)
permutest(beta)
TukeyHSD(beta)

boxplot(betadisper(wunifrac_dist,get_variable(ps.rarefied, "Phenotype")), las=2, main=paste0("Multivariate Dispersion Test uwunifrac "," pvalue = ", permutest(betadisper(uwunifrac_dist, get_variable(ps.rarefied, "Phenotype")))$tab$`Pr(>F)`[1]))
#uwunifrac_dist
anova(betadisper(uwunifrac_dist, rrsampledf1$Phenotype))
beta <- betadisper(uwunifrac_dist, rrsampledf1$Phenotype)
permutest(beta)
TukeyHSD(beta)

boxplot(betadisper(uwunifrac_dist,get_variable(ps.rarefied, "Phenotype")), las=2, main=paste0("Multivariate Dispersion Test uwunifrac "," pvalue = ", permutest(betadisper(uwunifrac_dist, get_variable(ps.rarefied, "Phenotype")))$tab$`Pr(>F)`[1]))
###save info
ggsave("C:/Users/57322/Desktop/cod/ord4_300.png", plot = ord4, width = 8, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("C:/Users/57322/Desktop/cod/ord3_300.png", plot = ord3, width = 8, height = 5, units = "in", dpi = 300, bg = "white")


arrange <- ggarrange(ay + 
                       theme(legend.position="none", strip.text = element_text(size = 16, colour = "black")) + font("legend.text", size = 16, color = "black", face = "bold")+ font("axis.title", size = 18, color = "black", face = "bold") + font("axis.text", size = 16, color = "black", face = "bold"), 
                     s + 
                       theme(strip.text = element_text(size = 12, colour = "black", face = "italic"))+ font("axis.title", size = 18, color = "black", face = "bold") + font("legend.text", color = "black", face = "italic", size = 14) + font("axis.text", size = 15, color = "black", face = "bold") + font("legend.title", size = 20, color = "black", face = "bold"),
                     ord3 + 
                       font("legend.text", size = 14, color = "black", face = "bold")+ font("axis.title", size = 18, color = "black", face = "bold")+ font("axis.text", size = 15, color = "black", face = "bold")+
                       theme(legend.position="none"), 
                     ord4 +
                       font("legend.title", size = 20, color = "black", face = "bold")+  font("legend.text", size = 14, color = "black", face = "bold")+ font("axis.title", size = 18, color = "black", face = "bold")+  font("axis.text", size = 15, color = "black", face = "bold") + font("subtitle", color = "black", face = "italic", size = 14))#, 

poiu <- as_ggplot(arrange) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C","D" ), size = 25,
                  x = c(0, 0, 0.42,0.42), y = c( 1, 0.48, 1, 0.48)) # Add labels
poiu
#save info
ggsave("C:/Users/57322/Desktop/cod/arrange_300.png", plot = poiu, width = 20, height = 12, units = "in", dpi = 300)

# Add labels
arranges <- ggarrange(ord3 + 
                        font("legend.text", size = 14, color = "black", face = "bold")+ font("axis.title", size = 18, color = "black", face = "bold")+ font("axis.text", size = 15, color = "black", face = "bold")+
                        theme(legend.position="none"), 
                      ord4 +
                        font("legend.title", size = 20, color = "black", face = "bold")+  font("legend.text", size = 14, color = "black", face = "bold")+ font("axis.title", size = 18, color = "black", face = "bold")+  font("axis.text", size = 15, color = "black", face = "bold") + font("subtitle", color = "black", face = "italic", size = 14))

#save info
ggsave("C:/Users/57322/Desktop/cod/arrangepot_300.png", plot = arranges, width = 20, height = 12, units = "in", dpi = 300)
#dev.off()
### ANOSIM

aswe<-anosim(wunifrac_dist,sample_data(ps.rarefied1)$Phenotype, permutations = 9999)
uni <-anosim(uwunifrac_dist,sample_data(ps.rarefied1)$Phenotype, permutations = 9999)


### ADONIS

# Adonis test
adonis(braycur_dist ~ Phenotype, data = rrsampledf1, permutations = 1000)
variable_group = get_variable(rrsampledf1, "Phenotype")
variable_group = anosim(uwunifrac_dist, variable_group)

uww <- adonis(uwunifrac_dist ~ sample_data(ps.rarefied1)$Phenotype, data = rrsampledf1, permutations = 9999)
variable_group = get_variable(rrsampledf1, "Phenotype")
variable_group = anosim(uwunifrac_dist, variable_group)

wun<-adonis(wunifrac_dist ~ sample_data(ps.rarefied1)$Phenotype, data = rrsampledf1, permutations = 9999)
variable_group = get_variable(rrsampledf1, "Phenotype")
variable_group = anosim(wunifrac_dist, variable_group)


### Counfounders
## Abundance Differential analysis DEseq2

###################################

physeqGenus = tax_glom(balea_2, "Genus")

physeqFamily = tax_glom(balea_2, "Family")

physeqClass = tax_glom(balea_2, "Class")
physeqPhylum = tax_glom(balea_2, "Phylum")
#########
ps_healthy <- ps_filter(physeqFamily, Phenotype != "Patients")
ps_Patients <- ps_filter(physeqFamily, Phenotype != "Healthy")
####################################################

heat.sample_H <- plot_taxa_heatmap(ps_healthy,
                                   subset.top = 20,
                                   VariableA = "Phenotype",
                                   heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdPu")))(100),
                                   transformation = "log10",
                                   show_rownames = TRUE,
                                   #legend = FALSE
)
#ggsave("C:/Users/57322/Desktop/cod/heat_healthy_fam.sample_300.png", plot = heat.sample_H[["plot"]], width = 20, height = 12, units = "in", dpi = 300)
heat.sample_P <- plot_taxa_heatmap(ps_Patients,
                                   subset.top = 20,
                                   VariableA = "Phenotype",
                                   heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdPu")))(100),
                                   transformation = "log10",
                                   show_rownames = TRUE,
                                   #legend = FALSE
)
heater_P <- plot_taxa_heatmap(physeqFamily,
                              subset.top = 20,
                              VariableA = "Phenotype",
                              heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdPu")))(100),
                              transformation = "log10",
                              show_rownames = TRUE,
                              #legend = FALSE
                              fontsize = 16
)

#plot(heat.sample)
ggsave("C:/Users/57322/Desktop/cod/heat_Patients_fam.sample_300.png", plot = heater_P[["plot"]], width = 20, height = 12, units = "in", dpi = 300)
###################




##############################################Deseq2 Abundance analysis

#taked/adapted from: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html and http://www.castrolab.org/isme/biodiversity/biodiversity.html#analisis-de-abundancias-y-visualizaciones


pipo<- phyloseq_to_deseq2(balea_2, ~ Phenotype)

# and running deseq standard analysis:
micro_16s_deseq <- DESeq(pipo)


# pulling out our results table, we specify the object, the p-value we are going to use to filter our results, and contrast (first naming the column, then the two groups)
micro_16s_deseq_difference <- results(micro_16s_deseq, alpha=0.01, contrast=c("Phenotype", "Patients", "Healthy"))

# summary command
summary(micro_16s_deseq_difference) 

# subset this table to only include significance level
sigtab_micro_16s_deseq_difference <- micro_16s_deseq_difference[which(micro_16s_deseq_difference$padj < 0.01), ]

# significantly differentially abundant
summary(sigtab_micro_16s_deseq_difference) 

# Pu together taxonomic annotations for a quick look 
sigtab_micro_16s_deseq_difference <- cbind(as(sigtab_micro_16s_deseq_difference, "data.frame"), as(tax_table(balea_2)[row.names(sigtab_micro_16s_deseq_difference), ], "matrix"))

# Sort that table by the baseMean column
yup<- sigtab_micro_16s_deseq_difference[order(sigtab_micro_16s_deseq_difference$baseMean, decreasing=T), ]
#plot(yup)
yup2<- sigtab_micro_16s_deseq_difference[order(sigtab_micro_16s_deseq_difference$baseMean), ]
# save info
write.table(format(yup2,scientific=FALSE, digits=4), file ="C:/Users/57322/Desktop/cod/deseq_zotus.csv", dec=",")

remove_these <- zotus_good_abu

#Now we find the indicies of the rows that need to be removed

rows_to_remove <- which(row.names(yup2) %in% remove_these)

#And use the same technique you were trying to use before to remove rows.

df <- yup2[rows_to_remove,]

##
#select the significative row from the table
eliminar <- c("Zotu35",  "Zotu66",  "Zotu31",  "Zotu59",  "Zotu75",  "Zotu11")
rows_to_remove <- which(row.names(df) %in% eliminar )
dfi <- df[rows_to_remove,]
###########Export info

sw<-stargazer(df,                 # Export txt
              summary = FALSE,
              type = "text",
              out = "C:/Users/57322/Desktop/cod/yup_solo_importantes.txt")
sw
#save info
pdf("C:/Users/57322/Desktop/cod/yup.pdf")       # Export PDF
grid.table(yup)
write.table(format(df,scientific=FALSE, digits=4), file ="C:/Users/57322/Desktop/cod/deseq_importantes.csv", dec=",")
#dev.off()
################### saving some other info
zotus_good_abu <-rownames(carbom_abund1@otu_table)
ash <- rownames(df)
remove_these2 <- ash

# rows that need to be removed
asii<- as.data.frame(objphil@otu_table)
rows_to_remove2 <- which(row.names(asii) %in% remove_these2)

#And use the same technique you were trying to use before to remove rows.

df2 <- asii[rows_to_remove2,]
sw2<-stargazer(df2,                 # Export txt
               summary = FALSE,
               type = "text",
               out = "C:/Users/57322/Desktop/cod/objphil_solo_importantes.txt")

################################ Plotting graphics

#####points
sigtabgen = subset(sigtab_micro_16s_deseq_difference, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Family = factor(as.character(sigtabgen$Family), levels=names(x))
# Genus order

deseq2plot<-ggplot(dfi, aes(y=Genus, x=log2FoldChange, color=Family)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

deseq2plot<-deseq2plot+ font("legend.title", size = 20, color = "black", face = "bold")+  font("legend.text", size = 14, color = "black", face = "italic")+ font("axis.title", size = 18, color = "black", face = "bold")+  font("axis.text", size = 15, color = "black", face = "italic") + font("subtitle", color = "black", face = "italic", size = 14) 
#save info
ggsave("C:/Users/57322/Desktop/cod/deseq2plot.png", plot =deseq2plot, width = 20, height = 12, units = "in", dpi = 300)
deseq2plot
########bars
deseq2plotbar<-ggplot(dfi, aes(y=Genus, x=log2FoldChange, color=Family)) + 
  geom_bar(aes(fill = Family),stat="identity", position=position_dodge())+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
deseq2plotbar<-deseq2plotbar+ font("legend.title", size = 20, color = "black", face = "bold")+  font("legend.text", size = 14, color = "black", face = "italic")+ font("axis.title", size = 18, color = "black", face = "bold")+  font("axis.text", size = 15, color = "black", face = "italic") + font("subtitle", color = "black", face = "italic", size = 14) 

deseq2plotbar
#save info
ggsave("C:/Users/57322/Desktop/cod/deseq2plotbar.png", plot =deseq2plotbar, width = 20, height = 12, units = "in", dpi = 300)

########### family ####

pipolin<- phyloseq_to_deseq2(physeqFamily, ~Phenotype)

# deseq standard analysis:
micro_16s_Fam_deseq <- DESeq(pipolin)

# pulling out our results table
micro_16s_Fam_deseq_difference <- results(micro_16s_Fam_deseq, alpha=0.01,  contrast=c("Phenotype", "Patients", "Healthy"))

# summary command
summary(micro_16s_Fam_deseq_difference) 
# specified significance level
sigtab_micro_16s_Fam_deseq_difference <- micro_16s_Fam_deseq_difference[which(micro_16s_Fam_deseq_difference$padj < 0.01), ]

# table only contains significantly differentially abundant
summary(sigtab_micro_16s_Fam_deseq_difference) 
# Put together
sigtab_micro_16s_Fam_deseq_difference <- cbind(as(sigtab_micro_16s_Fam_deseq_difference, "data.frame"), as(tax_table(physeqFamily)[row.names(sigtab_micro_16s_Fam_deseq_difference), ], "matrix"))

#  sort by the baseMean column
yup_fam<- sigtab_micro_16s_Fam_deseq_difference[order(sigtab_micro_16s_Fam_deseq_difference$baseMean, decreasing=T), ]
yup_fam
sw<-stargazer(yup_fam,                 # Export txt
              summary = FALSE,
              type = "text",
              out = "C:/Users/57322/Desktop/cod/familia_solo_importantes.txt")

write.table(format(yup_fam,scientific=FALSE, digits=4), file ="C:/Users/57322/Desktop/cod/familias.csv", dec=",")
###########para Phylum ####

pipolin_phy<- phyloseq_to_deseq2(physeqPhylum, ~Phenotype)


# running deseq standard analysis:
micro_16s_deseq <- DESeq(pipolin_phy)


# results table
micro_16s_deseq_difference <- results(micro_16s_deseq, alpha=0.01, contrast=c("Phenotype",  "Patients", "Healthy"))

# summary command
summary(micro_16s_deseq_difference) 

# table to only with significance level
sigtab_micro_16s_deseq_difference <- micro_16s_deseq_difference[which(micro_16s_deseq_difference$padj < 0.01), ]

# significantly differentially abundant
summary(sigtab_micro_16s_deseq_difference) 

# put together
sigtab_micro_16s_deseq_difference <- cbind(as(sigtab_micro_16s_deseq_difference, "data.frame"), as(tax_table(physeqPhylum)[row.names(sigtab_micro_16s_deseq_difference), ], "matrix"))

# sort table by the baseMean column
yup_phy<- sigtab_micro_16s_deseq_difference[order(sigtab_micro_16s_deseq_difference$baseMean, decreasing=T), ]
yup_phy

sw1<-stargazer(yup_phy,                 # Export txt
              summary = FALSE,
              type = "text",
              out = "C:/Users/57322/Desktop/cod/philum_solo_importantes.txt")
write.table(format(yup_phy,scientific=FALSE, digits=4), file ="C:/Users/57322/Desktop/cod/philum.csv", dec=",")

##################################################

######END