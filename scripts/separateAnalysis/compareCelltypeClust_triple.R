#The second setting of clustering will act as intermediates

library(Seurat)
library(networkD3)
library(tidyverse)

args  <- commandArgs(trailingOnly=T)
workdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/"
celltype <- args[1] #TCells, BCells, Monocytes, ImmuneCells, EpithelialCells, StromalCells
datatype <- args[2] #escc, gej, combined
se1_file <- paste0(workdir, args[3]) #rc_glmpca
se2_file <- paste0(workdir, args[4]) #sct_allhvg_hmy
se3_file <- paste0(workdir, args[5]) #sct_hmy
label1 <- args[6] #glmpca
label2 <- args[7] #hmy_allhvg
label3 <- args[8] #hmy

getClust <- function(se,label){
	cn <- colnames(se[[]])
	clust <- se[[cn[grep("_res.",cn)]]]
	colnames(clust)[1]=paste0(label,"_clust")
	clust[,1]=paste0(label,clust[,1])
	clust$cellid <- rownames(clust)
	if(grepl("G_",clust$cellid[1]) || grepl("E_",clust$cellid[1])){
		clust<-clust%>%separate(cellid,into=c("type","cellid"),sep="_",extra="merge")%>%select(-type)
	}
	return(clust)
}

se1 <- readRDS(file=se1_file)
clust1 <- getClust(se1,label1)
rm(se1)
se2 <- readRDS(file=se2_file)
clust2 <- getClust(se2,label2)
rm(se2)
se3 <- readRDS(file=se3_file)
clust3 <- getClust(se3,label3)
rm(se3)

merge1 <- merge(clust1,clust2)
merge1 <- as.data.frame(merge1%>%group_by_at(c(2,3))%>%summarise(count=n())%>%ungroup)
colnames(merge1)=c("source","target","count")
merge2 <- merge(clust2,clust3)
merge2 <- as.data.frame(merge2%>%group_by_at(c(2,3))%>%summarise(count=n())%>%ungroup)
colnames(merge2)=c("source","target","count")

common <- intersect(unique(merge1$target),unique(merge2$source))
merge1 <- merge1[merge1$target %in% common, ]
merge2 <- merge2[merge2$source %in% common, ]

##generate sankey plot
nodes <- data.frame(name=unique(c(
	as.character(merge1$source),
	as.character(merge1$target),
	as.character(merge2$target)
	)),stringsAsFactors = FALSE)

nodes$index <- 0:(nrow(nodes) - 1)

merge1$start <- as.integer(nodes$index[match(merge1$source, nodes$name)])
merge1$end <- as.integer(nodes$index[match(merge1$target, nodes$name)])
merge2$start <- as.integer(nodes$index[match(merge2$source, nodes$name)])
merge2$end <- as.integer(nodes$index[match(merge2$target, nodes$name)])
clust_merge <- rbind(merge1,merge2)

p <- sankeyNetwork(
  Links = clust_merge,
  Nodes = nodes,
  Source = "start",
  Target = "end",
  Value = "count",
  NodeID = "name",
  units = "cells",
  fontSize = 10,
  nodeWidth = 30,
  fontFamily = "Arial",
  iterations = 0,
  height = 200,
  width = 500
)
saveNetwork(p,paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/cellident/",
	celltype,"_",datatype,"_",label1,"_",label2,"_",label3,".sankey.html"),selfcontained=FALSE)
