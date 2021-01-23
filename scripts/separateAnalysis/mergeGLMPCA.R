suppressPackageStartupMessages({
        library(Seurat)
        library(future)
        library(future.apply)
	library(SeuratWrappers)
        library(glmpca)
})

args  <- commandArgs(trailingOnly=T)
celltype <- args[1] #TCells, BCells, Monocytes, ImmuneCells, EpithelialCells, StromalCells
nHVG <- as.numeric(args[2]) #2000 consistent with snakemake workflow GLMPCA
nworker <- min(as.numeric(args[3]),length(availableWorkers()))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
set.seed(1129)

outdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/Combined/"
gej_dir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustStab/"
escc_dir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/"

se <- readRDS(file=paste0(gej_dir,celltype,".rds"))
se1 <- readRDS(file=paste0(escc_dir,celltype,".rds"))

se <- merge(x=se, y=se1, add.cell.ids=c("G","E"), merge.data=FALSE)

ndims <- c(10,20)

m <- GetAssayData(se, slot = "counts", assay = "RNA")
devs <- scry::devianceFeatureSelection(m)
dev_ranked_genes <- rownames(se)[order(devs, decreasing = TRUE)]
for(n in ndims){
        #Sparse matrices are coerced to dense matrice for  minibatch='none'; If this exhausts memory, consider setting minibatch to 'stochastic' or 'memoized'
        se <- RunGLMPCA(se, features = head(dev_ranked_genes,n=nHVG), L = n, minibatch='stochastic',
                reduction.name=paste0('glmpca',n), reduction.key=paste0("GLMPC",n,"_"))
}

saveRDS(se, file=paste0(outdir, celltype, "_rc_glmpca.rds"))

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
