# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning

rm(list=ls()) #delete rm list of variables in environment given by ls()
setwd("/home/rj/ownCloud/DRIVE/CNV_Klebsielly/R")
# setwd("D:/ROBIN/CloudStation/CNV_Klebsielly/R")

library(gplots) #heatmap.2
library(ade4) #dist.binary
library(cluster) #daisy
library(RColorBrewer)

## Stromy2 CNVxMLST IMPORT
samples_table<-read.csv("tabulkaSAMPLES.csv",header = FALSE, sep = ";")
melts_table<-read.csv("tabulkaMLST.csv",header = FALSE, sep = ";")
samples_Labels <-read.csv("MELTtypyS.csv",header = FALSE, sep = ";",stringsAsFactors=FALSE)
samples_Labels <-unlist(samples_Labels) #convert list to character vector
rownames(samples_table) <-samples_Labels
rownames(melts_table) <-samples_Labels

## HCLUST SAMPLES
distance_matrix1<-dist.binary(samples_table, method = 5)
# distance_matrix1<-daisy(samples_table,metric="gower")
distance_matrix1[is.na(distance_matrix1)] <- 0
distance_matrix1[is.nan(distance_matrix1)] <- 0
sum(is.infinite(distance_matrix1))  # THIS SHOULD BE 0
distance_matrix1<- as.matrix(distance_matrix1)
heatmap.2(distance_matrix1,Rowv = TRUE,Colv = "Rowv",symm=TRUE)

##HCLUST MLST
# distance_matrix2<-dist.binary(melts_table, method = 10)
distance_matrix2<-daisy(melts_table,metric="gower")
distance_matrix2[is.na(distance_matrix2)] <- 0
distance_matrix2[is.nan(distance_matrix2)] <- 0
sum(is.infinite(distance_matrix2))  # THIS SHOULD BE 0
distance_matrix2<- as.matrix(distance_matrix2)
heatmap.2(distance_matrix2,Rowv = TRUE,Colv = TRUE,symkey=TRUE)


## MERGEDTRIANGLE MATRIX  
# distance_matrix1<-daisy(samples_table,metric="gower")
distance_matrix1<-dist.binary(samples_table, method = 5)

distance_matrix2<-daisy(melts_table,metric="gower")
matrix_combined <- distance_matrix1
matrix_combined[upper.tri(matrix_combined)] <- distance_matrix2[upper.tri(distance_matrix2)]

heatmap.2(as.matrix(matrix_combined),Rowv = TRUE,Colv = TRUE,symm = TRUE)

## https://bioinformatics.stackexchange.com/questions/4626/plotting-two-heatmaps-with-the-same-order-of-genes

#2
distance_matrix1<-daisy(samples_table,metric="gower")
distance_matrix2<-daisy(melts_table,metric="gower")

hr1 <-hclust(distance_matrix1, method = "average")
hr2 <-hclust(distance_matrix2, method = "average")

mergedcols <- merge(as.dendrogram(hr1),as.dendrogram(hr2))
mergeddata <- cbind(distance_matrix1, distance_matrix2)

# colors<- brewer.pal(n = 9, name = 'BuGn')
# myCol<- colorRampPalette(colors)(length(mergeddata))

hd<-dist(mergeddata, method = "euclidean")
# hd<-1-cor(t(mergeddata), method="pearson")
hd[is.na(hd)] <- 0
hd[is.nan(hd)] <- 0

hr <- hclust(as.dist(hd), method="average")

heatmap.2(mergeddata, Rowv=as.dendrogram(hr), Colv=as.dendrogram(mergedcols),
          distfun=function(x) dist(x, method="euclidean"),
          hclustfun=function(x) hclust(x, method="average")
)



# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# source("http://bioconductor.org/biocLite.R")
# biocLite("OmicsMarkeR")
# BiocManager::install(c("OmicsMarkeR"))
# library(BiocManager)
