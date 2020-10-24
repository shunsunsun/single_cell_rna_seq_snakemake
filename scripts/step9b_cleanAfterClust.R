#source /appsnew/source/R-4.0.2share.sh
#R.Version()
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scater))

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
picked_clust <- unname(unlist(strsplit(snakemake@params[["clust"]],':')))
randSeed <- 1129L
set.seed(randSeed)

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
se[['ident']] <- com
