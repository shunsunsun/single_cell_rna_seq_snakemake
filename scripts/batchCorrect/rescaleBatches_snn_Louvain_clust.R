#Ref: https://bioconductor.org/packages/devel/bioc/vignettes/batchelor/inst/doc/correction.html#5_Batch_rescaling
#need to run at source /appsnew/source/R.4.0.2share.sh
options(stringsAsFactors = F)
suppressPackageStartupMessages({
	library(Seurat)
	library(SingleCellExperiment)
	library(batchelor)
	library(future)
})

args  <- commandArgs(trailingOnly=T)
n_cpu <- as.numeric(args[1])

infile <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm.rds"
intermfile <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_rescaleRes.rds"
outfile <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_rescaleBatches.rds"
nworker <- min(n_cpu,length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
rseed=1129L
set.seed(rseed)

#inspired by SeuratWrapper::RunFastMNN
#batchelor::rescaleBatches returns a corrected matrix of per-gene log-expression values, wrapped in a SummarizedExperiment containin batch-related metadata. This function operates on a per-gene basis so there is no need to perform subsetting other than to improve speed (so the default value of features is NULL, or e.g., features = 2000)
RunRescaleBatches <- function(object.list, assay = NULL, features = NULL,verbose = TRUE, ...) {
  if (!all(sapply(X = object.list, FUN = inherits, what = 'Seurat'))) {
    stop("'object.list' must be a list of Seurat objects", call. = FALSE)
  }
  if (length(x = object.list) < 2) {
    stop("'object.list' must contain multiple Seurat objects for integration", call. = FALSE)
  }
  assay <- ifelse(is.null(assay), DefaultAssay(object = object.list[[1]]), assay)
  #message(paste0("Working on the ", assay, " assay"))
  for (i in 1:length(x = object.list)) {
    DefaultAssay(object = object.list[[i]]) <- assay
  }
  if (is.numeric(x = features)) {
    if (verbose) {
      message(paste("Computing", features, "integration features"))
    }
    features <- SelectIntegrationFeatures(
      object.list = object.list,
      nfeatures = features,
      assay = rep(assay,length(object.list))
    )
  }
  objects.sce <- lapply(
    X = object.list,
    FUN = function(x, f) {
      return(as.SingleCellExperiment(x = subset(x = x, features = f)))
    },
    f = features
  )
  integrated <- merge(
    x = object.list[[1]],
    y = object.list[2:length(x = object.list)]
  )
  out <- do.call(
    what = batchelor::rescaleBatches,
    args = c(
      objects.sce,
      list(...)
    )
  )
  saveRDS(out,file=intermfile)
  return(out)
}

if(!file.exists(outfile)){
    se <- readRDS(file=infile)
    if(!file.exists(intermfile)){
	    se@meta.data$orig.ident=sapply(strsplit(se@meta.data$orig.ident,"-"),"[",1)
	    corrected <- RunRescaleBatches(object.list = SplitObject(se,split.by = "orig.ident"))
    }else{
	    corrected <- readRDS(file=intermfile)
    }
    corrected <- Seurat::RunPCA(as.matrix(assay(corrected,"corrected")), verbose = FALSE) ##By default nPC=50
    saveRDS(corrected, file=gsub(".rds","_withPCA.rds",intermfile))
    #print(dim(se))
    #print(dim(corrected_se))
    se[['rspca']] <- corrected
    saveRDS(se, file=outfile)
}else{
    print(paste0("Present: ",outfile))
    se <- readRDS(file=outfile)
}
for (n in c(10,30,50)){
        for(k in c(20,50,80)){
                print(paste0("rescaledBatches pc dimension: ",n,"; snn neighborhood: ",k))
                se <- RunUMAP(se, reduction='rspca',dims=1:n,reduction.name=paste0("UMAP_rsbpca",n,"_snn",k), reduction.key = paste0("UMAPrsbpca",n,"snn",k,"_"),n.neighbors=k)
                se <- FindNeighbors(se,reduction='rspca',dims=1:n,k.param=k,graph.name=paste0("rsbpca",n,"_snn",k))
                #Enable method = "igraph" to avoid casting large data to a dense matrix
                se <- FindClusters(se, graph.name=paste0("rsbpca",n,"_snn",k),resolution=seq(0.8,1.2,by=0.2),random.seed=rseed, method="igraph")
        }
}
saveRDS(se, file=gsub(".rds","_clustered.rds",outfile))
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
