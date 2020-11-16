##--------Some commands to explore SingleR prediction result (denoted "pred")
#Note: (1) HSC -- haematopoietic stem cell; MEP -- megakaryocytic-erythroid progenitors giving rise to the cells that produce red blood cells and platelets
#(2) cell type labels before fine-tuning (first.labels), after fine-tuning (labels) and after pruning (pruned.labels)

colnames(pred)
head(sort(table(pred$labels), decreasing=TRUE))

#In pruned.labels, low-quality assignments are replaced with NA
#The default outlier-based pruning operation is ineffective if a particular cell label is constantly misassigned
to.remove <- is.na(pred$pruned.labels)
table(Label=pred$labels, Removed=to.remove)
sum(to.remove)

#For multi-reference classification, check the 'winning' reference
table(pred$reference)
#Check the result of individual reference
head(pred$orig.results$BPE$labels)
#Check the source reference for each label
table(Label=pred$labels, Reference=pred$reference)

#Subset cells
actual.hsc <- pred$labels[sce$protocol=="sorted hematopoietic stem cells" & sce$sample!="JC4"]
head(sort(table(actual.hsc), decreasing=TRUE))

#Diagnosis(may not be suitable for a large number of cells)
#plotScoreHeatmap(pred)

#compare cell identification and clustering results
#Adjusted rand index (ARI) >0.5 indicates reasonable consistency. However, higher resolution could lead to clusters nested within labels and therefore reduce ARI
library(bluster)
print(pairwiseRand(finalClust, pred$labels, mode="index")

library(Seurat)
library(dplyr)
se=readRDS(file="data/GEJ_EP_pipeCompflted_clustered_ann.rds")
colnames(se@meta.data)
df=se[[c("orig.ident","seurat_clusters","SingleR_label","SingleR_cluster_label")]]
table(df[df$SingleR_cluster_label %in% c("CD8+ Tcm","Monocytes","Plasma cells","MEP"),]$SingleR_label)%>%as.data.frame()%>%arrange(desc(Freq))%>%head
table(df[df$SingleR_cluster_label %in% c("CD8+ Tcm","Monocytes","Plasma cells"),]$orig.ident)%>%as.data.frame()%>%arrange(desc(Freq))%>%head
