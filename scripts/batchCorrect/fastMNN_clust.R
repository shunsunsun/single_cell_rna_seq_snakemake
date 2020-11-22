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

se <- readRDS(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchfastMNN.rds")
for (n in c(10,30,50)){
	for(k in c(20,50,80)){
		print(paste0("mnn dimension: ",n,"; snn neighborhood: ",k))
		se <- RunUMAP(se, reduction='mnn',dims=1:n,reduction.name=paste0("UMAP_mnn",n,"_snn",k), reduction.key = paste0("UMAPmnn",n,"snn",k,"_"),n.neighbors=k)
		se <- FindNeighbors(se,reduction='mnn',dims=1:n,k.param=k,graph.name=paste0("mnn",n,"_snn",k))
		#Enable method = "igraph" to avoid casting large data to a dense matrix
		se <- FindClusters(se, graph.name=paste0("mnn",n,"_snn",k),resolution=seq(0.8,1.2,by=0.2),random.seed=rseed, method="igraph")
	}
}
#objects with duplicate keys (offending key: UMAP_), setting key to 'umap_dim30_snn_'
#Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from UMAP_Dim10_SNN_ to UMAPDim10SNN_
#Warning: All keys should be one or more alphanumeric characters followed by an underscore '_', setting key to UMAPDim10SNN_
saveRDS(se, file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchfastMNN_clustered.rds")
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
