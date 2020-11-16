####---------Assign cell types based on scRNA-seq reference dataset: scHCL-------------
library(scHCL)
library(Seurat)
library(dplyr)

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
schclRes <- snakemake@params[[1]]

se <- readRDS(file=infile)
if(!file.exists(schclRes)){
	test=GetAssayData(object = se[["SCT"]], slot = "data")
	# scHCL has two parameters , single cell expression matrix(scdata) and the number of most similar cell types
	result <- scHCL(scdata = test, numbers_plot = 3)
	save(result, file=schclRes)
}else{
	load(file=schclRes)
}

#head(result$scHCL_probility)
df=result$scHCL_probility
df=df%>%arrange(Cell,desc(Score))%>%group_by(Cell)%>%top_n(n=1)%>%as.data.frame()
colnames(df)[2:3]=c("scHCL_label","scHCL_score")
df2=se@meta.data
df2$Cell=rownames(df2)
df2=merge(df2,df,by= "Cell")
se@meta.data <- df2
saveRDS(se, infile)

clust = df2%>%group_by(seurat_clusters)%>%summarise(clustSize=n())%>%as.data.frame()
scHCLpred = df2%>%group_by(seurat_clusters,scHCL_label)%>%count(scHCL_label)%>%arrange(seurat_clusters,desc(n))%>%group_by(seurat_clusters)%>%top_n(1)%>%as.data.frame()
singleRpred = df2%>%select(seurat_clusters,SingleR_cluster_label)%>%arrange(seurat_clusters)%>%as.data.frame()%>%unique()
compare = merge(clust,merge(scHCLpred,singleRpred))
colnames(compare)[4]="scHCLlabel_n"
saveRDS(compare,file=outfile)
