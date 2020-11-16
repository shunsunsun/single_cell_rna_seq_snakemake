library(dplyr)
df=readRDS(file="data/ESCC_EP_SingleR_scHCL_compare.rds")
df=df%>%mutate(scHCL_top_label_prop=scHCLlabel_n/clustSize)
df=df[,c(1,5,3,6)]
colnames(df)[c(1,3)]=c("clust","scHCL_cluster_top_label")
df$clust=as.numeric(df$clust)
df=df%>%arrange(clust)
write.table(df,col.names=T,row.names=F,sep='\t',quote=F)
