library(Seurat)
library(dplyr)
se=readRDS(file="data/ESCC_EP_pipeCompflted_clustered.rds")
#df=as.data.frame(table(se@meta.data$orig.ident))
#df%>%arrange(desc(Freq))
Idents(se) <- se$seurat_clusters

clusters=c(1,3:7,10,15)#GEJ_EP
clusters=c(1,3:8,11,12)#GEJ_IM
clusters=c(0,1:10,12,14) #ESCC_IM
clusters=c(3,5:8,11,14:16,18,25) #ESCC_EP

sub <- subset(se, idents=clusters) 
df2=as.data.frame(table(sub@meta.data$orig.ident))
df2=df2%>%arrange(desc(Freq))
tmp=table(sub@meta.data$seurat_clusters,sub@meta.data$orig.ident)

tmp1=tmp[as.character(clusters),]
df3=apply(tmp1,2,function(x){sum(x==0)})
samples=names(df3[df3==0])

tmp2=tmp[rownames(tmp)=="8",]
#apply(tmp2,2,function(x){sum(x==0)})
samples2=names(tmp2[tmp2>0])

colnames(df2)=c("SampleID","CellCount")
df2$Rank=rownames(df2)
df2=df2[df2$SampleID %in% intersect(samples,samples2), ]

tmp1=tmp1[,colnames(tmp1) %in% df2$SampleID]
