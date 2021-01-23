####---------Assign cell types based on scRNA-seq reference dataset: scHCL-------------
suppressPackageStartupMessages({
	library(scHCL)
	library(Seurat)
	library(dplyr)
})

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
tmpfile <- gsub(".rds","_scHCL.rda",infile)
singleRres <- gsub(".rds","_singleR_label.main.rda",infile)
outfile <- gsub(".rds","_autoAnnot.rds",infile)

se <- readRDS(file=infile)
if(!file.exists(tmpfile)){
	test=GetAssayData(object = se[["SCT"]], slot = "data")
	# scHCL has two parameters , single cell expression matrix(scdata) and the number of most similar cell types
	result <- scHCL(scdata = test, numbers_plot = 3)
	#save(result, file=tmpfile)
}else{
	load(file=tmpfile)
}

##Analyze the annotation and compare it with SingleR annotation
#head(result$scHCL_probility)
df=result$scHCL_probility
df=df%>%arrange(Cell,desc(Score))%>%group_by(Cell)%>%top_n(n=1)%>%as.data.frame()
colnames(df)[2:3]=c("scHCL_label","scHCL_score")
df2=se@meta.data
clustRes <- colnames(df2)[grepl("_res\\.",colnames(df2))]
df2 <- se[[clustRes]]
colnames(df2)="seurat_clusters"
df2$Cell=rownames(df2)
df2=merge(df2,df,by= "Cell")
rm(se,result)

clust = df2%>%group_by(seurat_clusters)%>%summarise(clustSize=n())%>%as.data.frame()
scHCLpred = df2%>%group_by(seurat_clusters,scHCL_label)%>%count(scHCL_label)%>%arrange(seurat_clusters,desc(n))%>%group_by(seurat_clusters)%>%top_n(1)%>%as.data.frame()

load(file=singleRres)
df2$SingleR_cluster_label <- clust.pred$labels[match(df2$seurat_clusters, rownames(clust.pred))]
singleRpred = df2%>%select(seurat_clusters,SingleR_cluster_label)%>%arrange(seurat_clusters)%>%as.data.frame()%>%unique()

compare = merge(clust,merge(scHCLpred,singleRpred))
colnames(compare)[4]="scHCLlabel_n"
saveRDS(compare,file=outfile)
