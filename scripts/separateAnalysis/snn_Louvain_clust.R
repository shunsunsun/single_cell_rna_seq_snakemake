args  <- commandArgs(trailingOnly=T)
infile <- args[1]
n_cpu <- as.numeric(args[2])
previous <- args[3] ##sct or cca or glmpca
rmTRBV <- args[4] #T or F
outPrefix <- args[5]

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
meta <- se@meta.data
previousRes <- colnames(meta)[grepl("res\\.", colnames(meta)) | grepl("seurat", colnames(meta))]
for(c in previousRes){
	se[[c]] <- NULL
}
print(infile)
if(previous=="sct"){
	se <- DietSeurat(se,scale.data = TRUE)
	DefaultAssay(se) <- "SCT"
	if(rmTRBV=="T"){
		hvg <- se@assays$SCT@var.features
		hvg <- setdiff(hvg, gene2rm)
		se <- RunPCA(object=se, features=hvg, verbose=F)  ##https://satijalab.org/seurat/v3.2/sctransform_vignette.html
	}else{
		se <- RunPCA(se, verbose = F) ##By default nPC=50
	}
}
if(previous=="cca"){
	DefaultAssay(se) <- "integrated"
	if(rmTRBV=="T"){
		hvg <- se@assays$integrated@var.features
		hvg <- setdiff(hvg, gene2rm)
		se <- ScaleData(se, features=hvg, verbose = F)
		se <- RunPCA(se, features=hvg, verbose = F)
	}else{
		se <- ScaleData(se, verbose = F)
		se <- RunPCA(se, verbose = F)
	}
}

if(!grepl("glmpca", previous)){
	tiff(file=paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/clustering/",
		outPrefix, "_pcHeatmap.tiff"), width=10, height=8, unit="in",res=100)
	DimHeatmap(se, dims=1:9, cells=500, balanced=T)
	dev.off()
	reduc="pca"
	ndims=c(20,40)
	label="pca"
}else{
	reduc=previous
	ndims=as.numeric(gsub("glmpca","",reduc))
	label="glmpca"
}

print(se[[reduc]],dim=1:9,nfeatures=20)
for (n in ndims){
	for(k in c(30,60,90)){
		print(paste0("pc dimension: ",n,"; snn neighborhood: ",k))
		#se <- RunUMAP(se, reduction=reduc,dims=1:n,reduction.name=paste0("UMAP_",label,n,"_snn",k), reduction.key = paste0("UMAP_",label,n,"snn",k,"_"),n.neighbors=k)
		se <- FindNeighbors(se,reduction=reduc,dims=1:n,k.param=k,graph.name=paste0(label,n,"_snn",k), verbose=F)
		#Enable method = "igraph" to avoid casting large data to a dense matrix
		se <- FindClusters(se,graph.name=paste0(label,n,"_snn",k),resolution=c(0.05,seq(0.1,0.9,by=0.2)),random.seed=rseed, method="igraph", verbose=F)
	}
}
if(rmTRBV=="T"){
	saveRDS(se, file=paste0(dirname(infile), "/", outPrefix, "_ignoreTRBV_clustered.rds"))
}else{
	saveRDS(se, file=paste0(dirname(infile), "/", outPrefix, "_clustered.rds"))
}
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
