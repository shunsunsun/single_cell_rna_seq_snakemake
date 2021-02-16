suppressPackageStartupMessages({
        library(Seurat)
        library(future)
        library(future.apply)
	library(SeuratWrappers)
        library(glmpca)
})

args  <- commandArgs(trailingOnly=T)
cancertype <- args[1] #ESCC or GEJ
nHVG <- as.numeric(args[2]) #2000 consistent with snakemake workflow GLMPCA
nworker <- min(as.numeric(args[3]),length(availableWorkers()))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
set.seed(1129)

workdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/FACS_based/"

se <- readRDS(file=paste0(workdir,cancertype,"_EP_merged.rds"))
se1 <- readRDS(file=paste0(workdir,cancertype,"_IM_merged.rds"))

se <- merge(x=se, y=se1, merge.data=FALSE)

ndims <- c(10,20)

m <- GetAssayData(se, slot = "counts", assay = "RNA")
devs <- scry::devianceFeatureSelection(m)
dev_ranked_genes <- rownames(se)[order(devs, decreasing = TRUE)]
for(n in ndims){
        #Sparse matrices are coerced to dense matrice for  minibatch='none'; If this exhausts memory, consider setting minibatch to 'stochastic' or 'memoized'
        se <- RunGLMPCA(se, features = head(dev_ranked_genes,n=nHVG), L = n, minibatch='stochastic',
                reduction.name=paste0('glmpca',n), reduction.key=paste0("GLMPC",n,"_"))
}

saveRDS(se, file=paste0(workdir, cancertype, "_rc_glmpca.rds"))

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
