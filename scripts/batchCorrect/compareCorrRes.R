args  <- commandArgs(trailingOnly=T)
n_cpu <- as.numeric(args[1])

suppressPackageStartupMessages({
	library(Seurat)
	library(future)
})
	
nworker <- min(n_cpu,length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
rseed=1129L
set.seed(rseed)

se <- readRDS(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm.rds")
se <- RunPCA(se, verbose = FALSE) ##By default nPC=50
for (n in seq(10,30,50)){
    print(n)
    se <- RunUMAP(se, reduction='pca',dims=1:n,reduction.name=paste0("UMAP_pca_dim",n), reduction.key = paste0("UMAPpcaDim",n,"_"))
    se <- FindNeighbors(se,reduction='pca',dims=1:n,graph.name=paste0("pca_dim",n))
    #Enable method = "igraph" to avoid casting large data to a dense matrix
    se <- FindClusters(se, graph.name=paste0("ccapca_dim",n),resolution=seq(0.4,1.2,by=0.2),random.seed=rseed, method="igraph")
}
saveRDS(se, file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustered.rds")
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
