# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning

rm(list=ls()) #delete rm list of variables in environment given by ls()
setwd("/home/rj/ownCloud/DRIVE/CNV_Klebsielly/R")

## CNV tree 1 IMPORT
data <-read.csv("tabulkaKP48.csv",header = TRUE, sep = ";")
CNVid <-read.csv("labels.csv",header = FALSE, sep = ";",stringsAsFactors=FALSE)
CNVid <-unlist(CNVid) #convert list to character vector

clusters <-c(unlist(data[2]),use.names=FALSE)
start<-c(unlist(data[5], use.names=FALSE))
stop<-c(unlist(data[6], use.names=FALSE))
coordinates <-cbind(start, stop) 
rownames(coordinates) <-CNVid

#HCLUST
distance_matrix<-dist(coordinates, method = "euclidean")
tree <-hclust(distance_matrix, method = "average")
# tree$labels

## 2 plot.dendrogram() function

plot(as.dendrogram(tree), type = "rectangle", 
     main = "CNV Coordinates Cluster Dendrogram", 
     ylab = "CNV", xlab = "Coordinates Euclidean Distance", horiz=FALSE)

plot(as.dendrogram(tree), ylim = c(1, 1e+02), xlim = c(1,62), horiz=FALSE) #zoom to axis limits

nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
plot(as.dendrogram(tree), type = "rectangle", 
     main = "CNV Coordinates Cluster Dendrogram", 
     ylab = "CNV", xlab = "Coordinates Euclidean Distance", 
     nodePar = nodePar, horiz = TRUE)

# 3 Phylogenetic trees
library("ape")

plot(as.phylo(tree), type = "phylogram", 
     main = "CNV Coordinates Cluster Dendrogram", 
     show.tip.label = TRUE,
     edge.color = "black", edge.width = 1, edge.lty = 1,
     tip.color = "black",
     cex = 0.6, label.offset = 0.5 )

# Fan
# Cut the dendrogram into clusters
colors = c("aquamarine", "violet", "greenyellow", "red",
             "cyan", "darkorange", "green", "royalblue",
             "magenta", "chartreuse", "purple", "pink",
             "deepskyblue", "darkviolet", "seagreen", "orangered")

png("rplot_CNV_fan.png", width = 1000, height = 1000, units = "px")

plot(as.phylo(tree), type = "fan", use.edge.length=FALSE,
       tip.color = colors[clusters],
       label.offset = 0.02, cex=1.75, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
dev.off()

##SILHOUETTE
library(cluster)

png("silhouette.png", width = 1000, height = 1000, units = "px")

sil<-silhouette(clusters, distance_matrix)
plot(sil)

dev.off()




## HEATMAP
library("RColorBrewer")
# display.brewer.all()  

distance_matrix<-as.matrix(distance_matrix)
# colors<- brewer.pal(6, "PuRd")#YlOrRd  BuGn PuRd
colors<- brewer.pal(n = 8, name = 'BuGn')
col<- colorRampPalette(colors)(length(distance_matrix))
# colorRampPalette(brewer.pal(8,"Dark2"))(length(distance_matrix))
heatmap(distance_matrix, col=col, symm=TRUE)  #heatmapa, need conversion to numerical matrix
