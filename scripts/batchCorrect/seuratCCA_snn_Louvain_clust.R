args  <- commandArgs(trailingOnly=T)
n_cpu <- as.numeric(args[1])
if(args[2]!=NULL){
	n_pc <- as.numeric(args[2])
}else{
	n_pc <- 50
}

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

se <- readRDS(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA.rds")
#Do NOT run the ScaleData function after integration
se <- RunPCA(se, verbose = FALSE, npcs=n_pc) ##By default npcs=50
for (n in c(10,30,50)){
	for(k in c(20,50,80)){
		print(paste0("pc dimension: ",n,"; snn neighborhood: ",k))
		se <- RunUMAP(se, reduction='pca',dims=1:n,reduction.name=paste0("UMAP_ccapca",n,"_snn",k), reduction.key = paste0("UMAPccapca",n,"snn",k,"_"),n.neighbors=k)
		se <- FindNeighbors(se,reduction='pca',dims=1:n,k.param=k,graph.name=paste0("ccapca",n,"_snn",k))
		#Enable method = "igraph" to avoid casting large data to a dense matrix
		se <- FindClusters(se, graph.name=paste0("ccapca",n,"_snn",k),resolution=seq(0.8,1.2,by=0.2),random.seed=rseed, method="igraph")
	}
}
#saveRDS(se, file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustered.rds")
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
