#Example: sbatch -p fat_icg -A gaog_g1 --qos=gaogfaticg -N 1 -n 1 -c 48 ./runRscript.sh integrateSubsets_CCA.R ImmuneCells 2000 cca 48

suppressPackageStartupMessages({
        library(Seurat)
        library(future)
        library(future.apply)
})

args  <- commandArgs(trailingOnly=T)
celltype <- args[1] #TCells, BCells, Monocytes, ImmuneCells, EpithelialCells, StromalCells
nFeature <- as.numeric(args[2])
reduc <- args[3] ##cca or rpca
nworker <- min(as.numeric(args[4]),length(availableWorkers()))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
set.seed(1129)

outdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/Combined/"
gej_dir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustStab/"
escc_dir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/"

obj.list <- list()
se <- readRDS(file=paste0(gej_dir,celltype,"_sctNorm_BatchHmy.rds"))
if(reduc=="cca"){
	se <- DietSeurat(se, assay=c("RNA","SCT"),scale.data=TRUE)
}else{
	#"rpca" requires precalcuated pca results
	se <- DietSeurat(se, assay=c("RNA","SCT"),dimreducs="pca",scale.data=TRUE)
}
obj.list[['gej']] <- se
se <- readRDS(file=paste0(escc_dir,celltype,"_sctNorm_BatchHmy.rds"))
if(reduc=="cca"){
	se <- DietSeurat(se, assay=c("RNA","SCT"),scale.data=TRUE)
}else{
	se <- DietSeurat(se, assay=c("RNA","SCT"),dimreducs="pca",scale.data=TRUE)
}
obj.list[['escc']] <- se
rm(se)

#Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = nFeature)

# Prepare the SCT list object for integration
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                anchor.features = integ_features, reduction = reduc, verbose = FALSE)

all_features <- Reduce(intersect, lapply(obj.list, rownames))
se <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", features.to.integrate = all_features, verbose=FALSE)

saveRDS(se, file=paste0(outdir, celltype, "_sctNorm_", toupper(reduc), ".rds"))

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
