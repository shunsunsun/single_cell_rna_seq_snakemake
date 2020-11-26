#need to run at source /appsnew/source/R.4.0.2share.sh
options(stringsAsFactors = F)
suppressPackageStartupMessages({
	library(sctransform)
	library(Seurat)
	library(sva)
	library(future)
})

args  <- commandArgs(trailingOnly=T)
n_cpu <- as.numeric(args[1])

infile <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm.rds"
outfile <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCombat_clustered.rds"
intermfile <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_combatCorrAssay.rds"
#covariates <- snakemake@params[["covariates"]] #concatenated by "+"
covariates <- NULL
nworker <- min(n_cpu,length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
rseed=1129L
set.seed(rseed)

###### break sparse matrix into trunk, convert into dense matrix and then concatenate
preprocess_big <- function(sparse_matrix,by=10000){
  n_col = ncol(sparse_matrix)
  n = floor(n_col/by)
  if(n<1){
    return(as.matrix(sparse_matrix))
  }else{
    for(i in 1:n){
      mat <- as.matrix(sparse_matrix[,((i-1)*by+1):(i*by)])
      if(i<2){
        res <- mat
      }else{
        res <- cbind(res,mat)
      }
      rm(mat)
    }
    if(n_col > n*by){
      res <- cbind(res,as.matrix(sparse_matrix[,(n*by+1):n_col]))
    }
  }
  rm(sparse_matrix)
  return(res)
}

if(!file.exists(outfile)){
    se <- readRDS(file=infile)
    df <- se@meta.data
    df$orig.ident=sapply(strsplit(df$orig.ident,"-"),"[",1)

    if(!file.exists(intermfile)){
	mod <- NULL
	if(!is.null(covariates)){
	    mod <- stats::model.matrix(stats::as.formula(paste0("~", covariates)), data = se@meta.data)
	}
	mat <- GetAssayData(object = se[["SCT"]], slot = "data") 

	flted <- preprocess_big(mat)
	print(dim(flted))
	row_max = matrixStats::rowMaxs(flted)
	row_min = matrixStats::rowMins(flted)
	if(sum(row_max==row_min)>0){
	  flted = flted[row_max != row_min,]
	}
	rv_genes <- which((apply(flted, 1, var)!=0) == "FALSE")
	rv_genes_names <- rownames(flted)[rv_genes]
	flted <- flted[!(rownames(flted) %in% rv_genes_names),]
	print(dim(flted))

	#source("/gpfs2/gaog_pkuhpc/users/liny/transcriptome_analysis/tcga_gtex/script/runCombat.r")
	#my_iter=15
	#res_assay <- ComBatWrapper(data.frame(flted), se@meta.data$orig.ident, mod, my_iter)
	#saveRDS(res_assay, file=outfile)

	print(dim(df))
	df <- df[colnames(flted), ]
	res_assay <- ComBat(dat=as.matrix(flted), batch = df$orig.ident, mod = mod, par.prior=TRUE, mean.only=FALSE, prior.plots=FALSE)
	saveRDS(res_assay, file=intermfile)
    }else{
	res_assay <- readRDS(file=intermfile)
    }
    se[['combatpca']] <- Seurat::RunPCA(res_assay, verbose = FALSE) ##By default nPC=50
    saveRDS(se, file=outfile)
}else{
    print(paste0("Present: ",outfile))
    se <- readRDS(file=outfile)
}
for (n in c(10,30,50)){
        for(k in c(20,50,80)){
                print(paste0("combat pc dimension: ",n,"; snn neighborhood: ",k))
                se <- RunUMAP(se, reduction='combatpca',dims=1:n,reduction.name=paste0("UMAP_combatpca",n,"_snn",k), reduction.key = paste0("UMAPcombatpca",n,"snn",k,"_"),n.neighbors=k)
                se <- FindNeighbors(se,reduction='combatpca',dims=1:n,k.param=k,graph.name=paste0("combatpca",n,"_snn",k))
                #Enable method = "igraph" to avoid casting large data to a dense matrix
                se <- FindClusters(se, graph.name=paste0("combatpca",n,"_snn",k),resolution=seq(0.8,1.2,by=0.2),random.seed=rseed, method="igraph")
        }
}
saveRDS(se, file=outfile)
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
