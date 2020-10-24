#source /appsnew/source/R-4.0.2share.sh
#R.Version()
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(bluster))

infile <- snakemake@input[[1]]
finiFlg <- snakemake@output[[1]]
outdir <- snakemake@params[[1]]
if(!dir.exists(outdir)){
        dir.create(outdir,recursive=T)
}
randSeed <- 1129L
set.seed(randSeed)

#Ref: https://cran.microsoft.com/web/packages/clustree/vignettes/clustree.html
drawClusTree <- function(df,outPrefix,reducs,param,w,h){
	for(rd in reducs){
	     if(!file.exists(paste0(outPrefix,"_",rd,".jpg"))){
		jpeg(paste0(outPrefix,"_",rd,".jpg"), width = w, height = h, units="cm", res=300)
		clus.tree.out <- clustree(df,prefix=paste0(rd,param),edge_width=0.8,node_alpha = 0.8) + 
			scale_color_brewer(palette = "Set1") + 
			theme(legend.position = "bottom") + 
			scale_edge_color_continuous(low = "grey80", high = "red")
		print(clus.tree.out)
		dev.off()
	     }	
	     #for(m in c("mitoCountRatio","cc_difference","log10_total_counts","log10_total_features","featcount_ratio")){
	     for(m in c("mitoCountRatio","cc_difference","featcount_ratio")){
		  if(!file.exists(paste0(outPrefix,"_",rd,"_",m,".jpg"))){
			jpeg(paste0(outPrefix,"_",rd,"_",m,".jpg"), width = w, height = h, units="cm", res=300)
			clus.tree.out <- clustree(df,prefix=paste0(rd,param),edge_width=0.8,node_alpha = 0.8,
				node_colour=m,node_colour_aggr="mean") + 
				theme(legend.position = "bottom") + 
				scale_edge_color_continuous(low = "grey80", high = "red")
			print(clus.tree.out)
                	dev.off()
		  }
             }
	}
}


#bluster::approxSilhouette
drawApproxSilhouette <- function(sce,outPrefix,param,reducs){
    cnames=colnames(colData(sce))
    for(rd in reducs){
	mat <- reducedDim(sce, toupper(rd))
	clustering <- cnames[grepl(paste0(rd,param),cnames)]
	plist=list()
	for(c in clustering){
	   sil.approx <- approxSilhouette(mat, clusters=sce[[c]])
           sil.data <- as.data.frame(sil.approx)
	   sil.data$closest <- factor(ifelse(sil.data$width > 0, sce[[c]], sil.data$other))
           sil.data$cluster <- factor(sce[[c]])
	   plist[[c]] <- ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
			ggbeeswarm::geom_quasirandom(method="smiley",size=.5) + 
			labs(x=paste0(c," cluster"),y="Approx. Silhouette width")
	}
	if(!file.exists(paste0(outPrefix,"_",rd,".jpg"))){
	    if(length(plist)>3){
		jpeg(paste0(outPrefix,"_",rd,".jpg"), width = 70, height = 40, units="cm", res=300)
		do.call("grid.arrange", c(plist, ncol=3))
		dev.off()
	    }else{
		jpeg(paste0(outPrefix,"_",rd,".jpg"), width = 55, height = 30, units="cm", res=300)
		do.call("grid.arrange", c(plist, ncol=2))
		dev.off()
	    }
	}
    }
}


if(grepl("seClust",infile)){
        se <- readRDS(file=infile)
	df=se@meta.data
        if(grepl("outly",infile)){
 	       df$log10_total_counts=log10(df$nCount_RNA)
               df$log10_total_features=log10(df$nFeature_RNA)
               df$featcount_ratio=df$log10GenesPerUMI
	       reducs <- 'pca'
        }else{
	       reducs <- c('pca','glmpca10','glmpca20')
	}
        drawClusTree(df,paste0(outdir,"/",gsub(".rds","ree",basename(infile))),reducs,"Clust_res.",35,20)
        drawApproxSilhouette(as.SingleCellExperiment(se),paste0(outdir,"/",gsub(".rds","_silhouette",basename(infile))),"Clust_res.",reducs)
        #DimPlot(se)
}else if(grepl("walkClust",infile)){
        load(file=infile)
	df=as.data.frame(colData(sce))
        if(grepl("outly",infile)){
               df$log10_total_counts=log10(df$nCount_RNA)
               df$log10_total_features=log10(df$nFeature_RNA)
               df$featcount_ratio=df$log10GenesPerUMI
	       reducs <- 'pca'
        }else{
	       reducs <- c('pca','glmpca10','glmpca20')
	}
        drawClusTree(df,paste0(outdir,"/",gsub(".rda","ree",basename(infile))),reducs,"_k.",35,10)
        drawApproxSilhouette(sce,paste0(outdir,"/",gsub(".rda","_silhouette",basename(infile))),"_k.",reducs)
        #plotReducedDim(sce, "TSNE", colour_by="label")
}
write(paste0("All done for ",infile), file=finiFlg)
