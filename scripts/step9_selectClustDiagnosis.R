#source /appsnew/source/R-4.0.2share.sh
#R.Version()
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MUDAN))

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
outdir <- snakemake@params[["outPlot"]]
if(!dir.exists(outdir)){
        dir.create(outdir,recursive=T)
}
picked_clust <- unname(unlist(strsplit(snakemake@params[["clust"]],':')))
randSeed <- 1129L
set.seed(randSeed)


# Visualization (including assessing the low-dimension space which the selected clustering was based on)
plotClust <- function(se, outPrefix, reducPrefix, res_colname){

	## Feature plot with QC metrics
	metrics=c("mitoCountRatio","cc_difference","featcount_ratio")
	p=FeaturePlot(se, reduction = paste0(reducPrefix,"_umap"), features = metrics, ncol=3, pt.size = .4, order = TRUE, min.cutoff = 'q10', label=T, label.size=3)
	if(!file.exists(paste0(outPrefix,"_qc.jpg"))){
		ggsave(file=paste0(outPrefix,"_qc.jpg"), plot = p, width = 32, height = 8, device="jpeg",units="cm",dpi=300)
	}

	## TSNE plot grouped and colored by clusters
	p1 <- DimPlot(se, reduction = paste0(reducPrefix, "_tsne"), group.by=res_colname, label=T) + NoLegend()
	#UMAP grouped and colored by clusters
	p2 <- DimPlot(se, reduction = paste0(reducPrefix, "_umap"), group.by=res_colname, label=T) + NoLegend()
	if(!file.exists(paste0(outPrefix,"_tsne_umap.jpg"))){
		ggsave(file=paste0(outPrefix,"_tsne_umap.jpg"), plot=p1+p2, width = 20, height=10, device="jpeg",units="cm",dpi=300)
	}

	## TSNE plot grouped and colored by cell cycle phases
	p3 <- DimPlot(se, reduction = paste0(reducPrefix, "_tsne"), group.by='Phase')
	if(!file.exists(paste0(outPrefix,"_tsne_clustVsCC.jpg"))){
		ggsave(file=paste0(outPrefix,"_tsne_clustVsCC.jpg"), plot = p1+p3, width = 20, height = 8, device="jpeg",units="in",dpi=300)
	}

	## UMAP grouped and colored by cell cycle phases
	p4 <- DimPlot(se, reduction = paste0(reducPrefix, "_umap"), group.by='Phase')
	if(!file.exists(paste0(outPrefix,"_umap_clustVsCC.jpg"))){
		ggsave(file=paste0(outPrefix,"_umap_clustVsCC.jpg"), plot = p2+p4, width = 20, height = 8, device="jpeg",units="cm",dpi=300)
	}
	
	## TSNE and UMAP split and colored by cell cycle phases
	p5 <- DimPlot(se, reduction = paste0(reducPrefix, "_tsne"), group.by='Phase', split.by='Phase')
	p6 <- DimPlot(se, reduction = paste0(reducPrefix, "_umap"), group.by='Phase', split.by='Phase')
	if(!file.exists(paste0(outPrefix,"_cc_tsne_umap.jpg"))){
		ggsave(file=paste0(outPrefix,"_cc_tsne_umap.jpg"), plot = p5/p6, width = 30, height = 20, device="jpeg",units="cm",dpi=300)
	}

	#save(p,p1,p2,p3,p4,p5,p6,file=paste0(outPrefix,"_ggplot2Objs.rda"))
}

checkClustStability <- function(cd, com, matnorm, outfile, verbose=FALSE){
	##MUDAN::getStableClusters
	stable <- getStableClusters(cd, com, matnorm, verbose=verbose, plot=FALSE)
	com.sub <- na.omit(stable$com)
	com=as.numeric(unique(com))
	com.sub=as.numeric(unique(com.sub))
	write(paste0("There are ", length(com), " clusters: ", paste0(com[order(com)], collapse=',')), file=outfile)
	write(paste0("Among them, ", length(com.sub), " stable ones: ", paste0(com.sub[order(com.sub)], collapse=',')), file=outfile, append=T)
}


##Recall that dimension reduction was conducted on a Seurat object, which was then converted into a SingleCellExperiment object for walk trap clustering, 
##So we only need to add the walk trap clustering result into the meta.data of the Seurat object and then would be able to use the visusalization tools of Seurat
if(picked_clust[1]=="seurat"){
        se <- readRDS(file=infile)
	clustRes=se[[paste0(picked_clust[2],"Clust_res.",picked_clust[3])]]
	com=clustRes[,1]
	names(com)=rownames(clustRes)
	com=factor(com)
	#finalClustLabel <- c(paste0(picked_clust[1],'_K50'), picked_clust[2], paste0('res',picked_clust[3]))
}else if(picked_clust[1]=="walktrap"){
	se <- readRDS(file=infile)
        load(file=gsub("se","walk",gsub("rds","rda",infile)))
	df=colData(sce)
	com=df[, paste0(picked_clust[2],"_k.",picked_clust[3])]
	names(com)=rownames(df)
	com=factor(com)
	#finalClustLabel <- c(picked_clust[1], picked_clust[2], paste0('K',picked_clust[3]))
}
#se <- AddMetaData(se, metadata = com, col.name=paste(finalClustLabel, collapse='_'))
Idents(object = se) <- com
se[['seurat_clusters']] <- com
se[['ident']] <- com
#print(paste(unique(se[['ident']]),collapse=','))
#plotClust(se,paste0(outdir,"/",gsub("seClust.rds","finalClust",basename(infile))),picked_clust[2],'ident')
cd=GetAssayData(object = se[["RNA"]], slot = "counts")
matnorm=GetAssayData(object = se[["SCT"]], slot = "data")
if(!file.exists(outfile)){
	checkClustStability(cd,com,matnorm,outfile,TRUE)
}
