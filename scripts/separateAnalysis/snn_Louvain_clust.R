args  <- commandArgs(trailingOnly=T)
infile <- args[1]
n_cpu <- as.numeric(args[2])
previous <- args[3] ##sct or cca
rmTRBV <- args[4] #T or F

#gene2rm <- c("TRBV20-1","TRBV11-2","TRBV7-2")
gene2rm <- "TRBV20-1"

suppressPackageStartupMessages({
	library(Seurat)
	library(future)
})
	
nworker <- min(n_cpu,length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 70*1024^3)
rseed=1129L
set.seed(rseed)

#se <- readRDS(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm.rds")
se <- readRDS(file=infile)
if(previous=="sct"){
	se <- DietSeurat(se,scale.data = TRUE)
	DefaultAssay(se) <- "SCT"
	if(rmTRBV=="T"){
		hvg <- se@assays$SCT@var.features
		hvg <- setdiff(hvg, gene2rm)
		se <- RunPCA(object=se, features=hvg, verbose=FALSE)  ##https://satijalab.org/seurat/v3.2/sctransform_vignette.html
	}else{
		se <- RunPCA(se, verbose = FALSE) ##By default nPC=50
	}
}
if(previous=="cca"){
	DefaultAssay(se) <- "integrated"
	if(rmTRBV=="T"){
		hvg <- se@assays$integrated@var.features
		hvg <- setdiff(hvg, gene2rm)
		se <- ScaleData(se, features=hvg, verbose = FALSE)
		se <- RunPCA(se, features=hvg, npcs = 30, verbose = FALSE)
	}else{
		se <- ScaleData(se, verbose = FALSE)
		se <- RunPCA(se, npcs = 30, verbose = FALSE)
	}
}
for (n in c(10,30)){
	for(k in c(20,50,80)){
		print(paste0("pc dimension: ",n,"; snn neighborhood: ",k))
		se <- RunUMAP(se, reduction='pca',dims=1:n,reduction.name=paste0("UMAP_pca",n,"_snn",k), reduction.key = paste0("UMAPpca",n,"snn",k,"_"),n.neighbors=k)
		se <- FindNeighbors(se,reduction='pca',dims=1:n,k.param=k,graph.name=paste0("pca",n,"_snn",k))
		#Enable method = "igraph" to avoid casting large data to a dense matrix
		se <- FindClusters(se, graph.name=paste0("pca",n,"_snn",k),resolution=seq(0.2,1.2,by=0.2),random.seed=rseed, method="igraph")
	}
}
if(rmTRBV=="T"){
	saveRDS(se, file=gsub(".rds", "_ignoreTRBV_clustered.rds", infile))
}else{
	saveRDS(se, file=gsub(".rds", "_clustered.rds", infile))
}
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
