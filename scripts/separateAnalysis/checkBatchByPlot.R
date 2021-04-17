suppressPackageStartupMessages({
	library(Seurat)
	library(ggplot2)
	library(patchwork)
	library(future)
	library(tidyr)
})

randseed=1129L

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
cohort <- args[2] #e.g., ESCCjiaoda_TCells_sct_refCCA
plotfile <- paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/batchCorrect/", cohort, "_dimplot_batchcheck.pdf")
nworker <- min(as.numeric(args[3]),length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 20*1024^3)
reduc <- args[4] #pca or glmpca10 or glmpca20

se <- readRDS(file=infile)
if(grepl("jiaoda", cohort)){
	print("Jiaoda data")
	se$orig.ident <- gsub("_filtered_feature_bc_matrix","",se$orig.ident)
	se$type <- ifelse(grepl("A$", se$orig.ident), "T", "N")
	se$patient <- sapply(se$orig.ident, function(x){a=unlist(strsplit(x,"_"))[2];return(a)})
	se$patient <- gsub("A","",gsub("B","",se$patient))
}else if(grepl("TDN", cohort)){ #TDN data
	print("ESCC_TDN")
	meta <- se@meta.data
	meta$cid <- rownames(meta)
	clinic <- read.table(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_14patients_TDN_clustStab/sample.cellno.stage.txt",header=T)
	clinic <- clinic[,c(1,2)]
	colnames(clinic)[1:2] <- c("orig.ident","p_sample")
	clinic <- clinic %>% separate(p_sample, c("patient","type"))
	meta=merge(meta,clinic)
	rownames(meta)=meta$cid
	se@meta.data <- meta
	se$cid <- NULL	
}else{
	print("first sc cohort")
	se$type <- ifelse(grepl("N", se$orig.ident), "N", "T")
	se$patient=gsub("T","",gsub("N","",se$orig.ident))
}
if(reduc=="pca"){
	if(!grepl("_sct", cohort)){
		se <- ScaleData(se, verbose = FALSE)
	}
	# Run PCA
	se <- RunPCA(object = se)
	# Run UMAP
	se <- RunUMAP(se, dims=1:40, reduction = "pca")
}else{
	n <- as.numeric(gsub("glmpca","",reduc))
	se <- RunUMAP(se, dims=1:n, reduction = reduc)
}

p1 = DimPlot(se, group.by="orig.ident", shuffle=T, seed=randseed, label=F, reduction="umap")
p2 = DimPlot(se, group.by="patient", shuffle=T, seed=randseed, label=F, reduction="umap")
p3 = DimPlot(se, group.by="type", shuffle=T, seed=randseed, label=F, reduction="umap")
p4 = plot_spacer()
p = (p1 + p4) / (p2 + p3)
p = p + plot_annotation(tag_levels = 'A') 
ggsave(filename=plotfile, plot=p, width = 10, height = 8, units="in",dpi=300)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
