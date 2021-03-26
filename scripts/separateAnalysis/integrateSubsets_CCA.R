#Example: sbatch -p fat_icg -A gaog_g1 --qos=gaogfaticg -N 1 -n 1 -c 48 ./runRscript.sh integrateSubsets_CCA.R 

suppressPackageStartupMessages({
        library(Seurat)
        library(future)
        library(future.apply)
})

args  <- commandArgs(trailingOnly=T)
obj1_file <- args[1]
obj2_file <- args[2]
resfile <- args[3]
nFeature <- as.numeric(args[4]) #2000
reduc <- args[5] ##cca or rpca
nworker <- min(as.numeric(args[6]),length(availableWorkers()))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
set.seed(1129L)

obj.list <- list()
se <- readRDS(file=obj1_file)
se <- DietSeurat(se, assay=c("RNA","SCT"),scale.data=TRUE)
obj.list[[1]] <- se
rm(se)
se <- readRDS(file=obj2_file)
se <- DietSeurat(se, assay=c("RNA","SCT"),scale.data=TRUE)
obj.list[[2]] <- se
rm(se)

#Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = nFeature)

# Prepare the SCT list object for integration
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                anchor.features = integ_features, reduction = reduc, verbose = FALSE)

all_features <- Reduce(intersect, lapply(obj.list, rownames))
se <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", features.to.integrate = all_features, verbose=FALSE)

if(!dir.exists(dirname(resfile))){
        dir.create(dirname(resfile),recursive=T)
}
saveRDS(se, file=resfile)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
