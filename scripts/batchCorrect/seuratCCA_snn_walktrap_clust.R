args  <- commandArgs(trailingOnly=T)
n_cpu <- as.numeric(args[1])

suppressPackageStartupMessages({
	library(Seurat)
        library(scater)
        library(scran)
        library(BiocParallel)
})
	
nworker <- min(n_cpu,parallel::detectCores()-1)
cat(sprintf("Use %d workders\n",nworker))
rseed=1129L
set.seed(rseed)

se <- readRDS(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustered.rds")
se <- DietSeurat(se, assays='integrated', dimreducs='pca')
sce <- as.SingleCellExperiment(se)
for(k in c(20,50,80)){
	print(paste0("snn neighborhood: ",k))
	g <- buildSNNGraph(sce, type="jaccard", k=k, use.dimred = 'PCA',BPPARAM=MulticoreParam(nworker))
        sce[[paste0("snn",k)]] <- igraph::cluster_walktrap(g)$membership
}
saveRDS(sce, file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_walktrap_clustered.rds")
