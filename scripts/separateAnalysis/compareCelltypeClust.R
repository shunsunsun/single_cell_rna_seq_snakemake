library(Seurat)
library(networkD3)
library(tidyverse)

args  <- commandArgs(trailingOnly=T)
workdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/"
celltype <- args[1] #TCells, BCells, Monocytes, ImmuneCells, EpithelialCells, StromalCells
datatype <- args[2] #escc, gej, combined
se1_file <- paste0(workdir, args[3])
se2_file <- paste0(workdir, args[4])
label1 <- args[5] ##label example: glmpca, harmony
label2 <- args[6]

se1 <- readRDS(file=se1_file)
se2 <- readRDS(file=se2_file)

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
clust2 <- getClust(se2,label2)
clust1 <- getClust(se1,label1)
clust_merge <- merge(clust1,clust2)
clust_merge <- as.data.frame(clust_merge%>%group_by_at(c(2,3))%>%summarise(count=n())%>%ungroup)

##generate sankey plot
nodes <- data.frame(name=unique(c(
	as.character(clust_merge[,1]),
	as.character(clust_merge[,2])
	)),stringsAsFactors = FALSE)

nodes$index <- 0:(nrow(nodes) - 1)

clust_merge$source <- as.integer(nodes$index[match(clust_merge[,1], nodes$name)])
clust_merge$target <- as.integer(nodes$index[match(clust_merge[,2], nodes$name)])

p <- sankeyNetwork(
  Links = clust_merge,
  Nodes = nodes,
  Source = "target",
  Target = "source",
  Value = "count",
  NodeID = "name",
  units = "cells",
  fontSize = 10,
  nodeWidth = 30,
  fontFamily = "Arial",
  iterations = 0,
  height = 200,
  width = 300
)
saveNetwork(p,paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/cellident/",
	celltype,"_",datatype,"_",label1,"_",label2,".sankey.html"),selfcontained=FALSE)
