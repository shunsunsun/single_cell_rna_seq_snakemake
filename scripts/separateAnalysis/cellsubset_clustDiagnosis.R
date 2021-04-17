suppressPackageStartupMessages({
        library(Seurat)
	library(clustree)
	library(ggplot2)
	library(patchwork)
})

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
cohort <- args[2]
step <- as.numeric(args[3])
outdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/clustering/"
#se <- readRDS(file="../../data/ESCC_QCed_sctNorm_BatchHmy_clustStab/TCells_rc_glmpca20_k100_res0.15_C01_withSCT_noBatchCorr_clustered.rds")
se <- readRDS(file=infile)

if(step==1){
  for(pc in c(10,30)){
     for(k in c(30,60,90)){
	label <- paste0("pca",pc,"_snn",k)
	p <- clustree(se, prefix=paste0(label,"_res."),edge_width=0.8,node_alpha = 0.8) + 
		scale_color_brewer(palette = "Set1") + theme(legend.position = "bottom") + 
		scale_edge_color_continuous(low = "grey80", high = "red")
	ggsave(plot=p, filename=paste0(outdir, cohort, "_", gsub("clustered.rds",paste0(label,"_clustree.pdf"),basename(infile))), 
		width=10, height=10, units="in", dpi=300)
     }
  }
}

if(step==2){
  nPC <- args[4]
  k <- args[5]
  res <- args[6]
  meta <- se@meta.data
  cols <- colnames(meta)[grepl("_res\\.", colnames(meta)) & !(grepl(paste0("pca",nPC,"_snn",k,"_res.",res),colnames(meta)))]
  for(n in cols){
	se[[n]] <- NULL
  }
  se$seurat_clusters <- NULL
  reduc <- names(se)[grepl("UMAP_pca", names(se)) & !grepl(paste0("UMAP_pca", nPC, "_snn", k), names(se))]
  for(n in reduc){
	se[[n]] <- NULL
  }
  graph <- names(se)[grepl("pca", names(se)) & !grepl(paste0("pca",nPC,"_snn",k), names(se))]
  for(g in graph){
	se[[g]] <- NULL
  }
  saveRDS(se, file=infile)
}

if(step==3){
	randseed=1129L
	reduc <- names(se)[grepl("UMAP_pca",names(se))]
        group <- colnames(se[[]])[grepl("_res\\.",colnames(se[[]]))]
	p1 <- DimPlot(se, reduction=reduc, group.by=group, label=T) + NoLegend()
	p2 <- DimPlot(se, reduction=reduc, group.by="orig.ident", shuffle=T, seed=randseed, label=F) + NoLegend()
	#p3 <- DimPlot()
	p=p1+p2
	ggsave(plot=p,filename=paste0(outdir,cohort, "_",gsub("noBatchCorr_clustered.rds","dimplot.pdf",basename(infile))),
		 width=10, height=4, units="in", dpi=300)
}

