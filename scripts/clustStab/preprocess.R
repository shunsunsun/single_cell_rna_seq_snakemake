suppressPackageStartupMessages({
	library(Seurat)
	library(tidyverse)
	library(future)
})

## see https://bitbucket.org/snakemake/snakemake/issues/917/enable-stdout-and-stderr-redirection
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

seurat_obj<- readRDS(snakemake@input[[1]])
k <- snakemake@wildcards[["k"]]
resolution <- snakemake@params[["resolution"]]
pc.use<- snakemake@wildcards[["pc"]]
nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
output <- snakemake@output[[1]]

cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
rseed=1129L
set.seed(rseed)

PreprocessSubsetData<- function(object,
                                variable.features.n = 3000,
                                num.pc = 50,
                                pc.use = NULL,
                                workers = 2,
                                score.thresh = 1e-5,
                                sig.pc.thresh = 0.05,
                                n.start = 100,
                                nn.eps = 0,
                                resolution = seq(0.8,1.2,by=0.2),
                                k.param = 30,
                                ...){



	object <- readRDS(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA.rds")
	#Do NOT run the ScaleData function after integration
	object <- RunPCA(object, verbose = FALSE, npcs=num.pc) ##the default value of npcs = the default value of num.pc = 50
	
        if (is.null(pc.use)){
                object<- JackStraw( object = object, num.replicate = 100, dims = num.pc)

                object <- ScoreJackStraw(object = object, dims = 1:num.pc, score.thresh = score.thresh)

                PC_pvalues<- object@reductions$pca@jackstraw@overall.p.values

                ## determin how many PCs to use.
                pc.use<- min(which(PC_pvalues[,"Score"] > sig.pc.thresh)) -1

        }
	
	# add significant pc number to metadata, need to have names same as the cells
        pc.use.meta<- rep(pc.use, length(colnames(object)))
        names(pc.use.meta)<- colnames(object)
        object<- AddMetaData(object = object, metadata = pc.use.meta, col.name = "pc.use")
	
	object <- FindNeighbors(object,reduction='pca',dims=1:pc.use,k.param=k.param,verbose=FALSE,force.recalc=TRUE)
	
	#Enable method = "igraph" to avoid casting large data to a dense matrix
	object <- FindClusters(object, reduction.type='pca', n.start=n.start, resolution=resolution,
		random.seed=rseed, method="igraph",verbose=FALSE)
	
	return(object)
}


PreprocessSubsetData_pars<- snakemake@params[["PreprocessSubsetData_pars"]]
## this is not subsetted data, but the PreprocessSubsetData function can be used as well for any seurat object
seurat_obj<- eval(parse(text=paste("PreprocessSubsetData", "(", "seurat_obj,", "k.param=", k, ",",
	"pc.use=", pc.use, ",", PreprocessSubsetData_pars, ")")))
#saveRDS(seurat_obj, file = paste0("full_sample_preprocess/full_sample_", "k_", k, "_PC_", pc.use, ".rds"))
saveRDS(seurat_obj,file=output)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
