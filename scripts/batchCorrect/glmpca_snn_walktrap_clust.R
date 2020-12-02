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

se <- readRDS(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_glmpcaReduct_clustered.rds")
se <- DietSeurat(se, assays='RNA', dimreducs=c('glmpca10','glmpca20'))
sce <- as.SingleCellExperiment(se)

for (n in c(10,20)){
	for(k in c(20,50,80)){
		print(paste0("glmpca dimension: ",n,"; snn neighborhood: ",k))
		g <- buildSNNGraph(sce, type="jaccard", k=k, use.dimred = paste0('GLMPCA',n),BPPARAM=MulticoreParam(nworker))
	        sce[[paste0("glmpca",n,"_snn",k)]] <- igraph::cluster_walktrap(g)$membership
	}
}
saveRDS(sce, file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_glmpcaReduct_walktrap_clustered.rds")
