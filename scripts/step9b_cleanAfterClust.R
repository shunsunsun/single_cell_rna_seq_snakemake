#source /appsnew/source/R-4.0.2share.sh
#R.Version()
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scater))

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
picked_clust <- unname(unlist(strsplit(snakemake@params[["clust"]],'_')))

##Recall that dimension reduction was conducted on a Seurat object, which was then converted into a SingleCellExperiment object for walk trap clustering, 
##So we only need to add the walk trap clustering result into the meta.data of the Seurat object and then would be able to use the visusalization tools of Seurat
se <- readRDS(file=infile)
if(picked_clust[1]=="seurat"){
	clustRes=se[[paste0(picked_clust[2],"Clust_res.",picked_clust[3])]]
	com=clustRes[,1]
	names(com)=rownames(clustRes)
	com=factor(com)
}else if(picked_clust[1]=="walktrap"){
        load(file=gsub("se","walk",gsub("rds","rda",infile)))
	df=colData(sce)
	com=df[, paste0(picked_clust[2],"_k.",picked_clust[3])]
	names(com)=rownames(df)
	com=factor(com)
}
Idents(object = se) <- com
se[['seurat_clusters']] <- com
se <- DietSeurat(se, counts = TRUE, data = TRUE, scale.data = TRUE,graph = paste0(picked_clust[2],"Clust"),
        dimreducs = c(picked_clust[2],paste0(picked_clust[2],'_tsne'),paste0(picked_clust[2],'_umap')))
col2rm <- colnames(se[[]])[grepl("_res.",colnames(se[[]]))] 
##https://github.com/satijalab/seurat/issues/2017#issuecomment-524418929
for(c in col2rm){
	se[[c]] <- NULL
}   
saveRDS(se, file=outfile)
