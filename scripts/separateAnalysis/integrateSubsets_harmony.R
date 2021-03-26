#Example: sbatch -p cn_icg -A gaog_g1 --qos=gaogcnicg -N 1 -n 1 -c 30 ./runRscript.sh integrateSubsets_harmony.R BCells IM T 30

suppressPackageStartupMessages({
        library(Seurat)
	library(harmony)
        library(future)
        library(future.apply)
})

args  <- commandArgs(trailingOnly=T)
celltype <- args[1] #TCells, BCells, Monocytes, ImmuneCells, EpithelialCells, StromalCells
facs <- args[2] #IM or EP
regressNum <- args[3]
nworker <- min(as.numeric(args[4]),length(availableWorkers()))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
set.seed(1129L)

outdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/Combined/"
gej_dir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustStab/"
escc_dir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/"

escc_facs <- readRDS(file = paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/FACS_based/ESCC_",toupper(facs),"_pipeCompflted_clustered_SingleRann.rds"))
gej_facs <- readRDS(file = paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/FACS_based/GEJ_",toupper(facs),"_pipeCompflted_clustered_SingleRann.rds"))


filterCells <- function(obj, facs_data){
	cells.used <- unname(facs_data[["Cell", drop=TRUE]])
	obj <-subset(obj, cells=cells.used)		
}

se <- readRDS(file=paste0(gej_dir,celltype,".rds"))
se <- filterCells(se, gej_facs)
se1 <- readRDS(file=paste0(escc_dir,celltype,".rds"))
se1 <- filterCells(se1, escc_facs)

se <- merge(x=se, y=se1, merge.data=FALSE)
rm(se1)

if(regressNum=="T"){
    regVars=c('mitoCountRatio', 'nFeature_RNA', 'nCount_RNA')
}else{
    regVars=c('mitoCountRatio')
}

se <- SCTransform(se, variable.features.n=3000, vars.to.regress = c(regVars, 'cc_difference'))
se <- RunPCA(object=se, verbose=T)  ##https://satijalab.org/seurat/v3.2/sctransform_vignette.html
se@meta.data$orig.ident=sapply(strsplit(se@meta.data$orig.ident,"-"),"[",1)
se <- RunHarmony(object=se, assay.use = "SCT", reduction = "pca", dims.use = 1:50, group.by.vars = "orig.ident", plot_convergence = FALSE)

saveRDS(se, file=paste0(outdir, celltype, "_sctNorm_BatchHmy.rds"))

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
