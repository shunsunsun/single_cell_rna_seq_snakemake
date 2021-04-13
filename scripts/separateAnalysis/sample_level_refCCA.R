options(stringsAsFactors = F)
suppressPackageStartupMessages({
        library(sctransform)
        library(Seurat)
	library(dplyr)
        library(future)
        library(future.apply)
})

refs <- c("GSM4317409_S133A_filtered_feature_bc_matrix","GSM4317412_S149A_filtered_feature_bc_matrix","GSM4317413_S150A_filtered_feature_bc_matrix") #For ESCC_jiaoda
#refs <- c("LZE8T","LZE7T","LZE6N") #For ESCC_14TDN

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
outfile <- args[2]
targetAssay <- "SCT"
nFeature <- 2000
mode <- "cca"

# prevent small datasets from not having enough neighbors (k) to use when filtering anchors
# https://github.com/satijalab/seurat/issues/997#issuecomment-505207415
k_filter <- 200

nworker <- min(as.numeric(args[3]),length(availableWorkers()))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
set.seed(1129L)

se <- readRDS(file=infile)
DefaultAssay(se) <- "SCT"
se <- DietSeurat(se, scale.data=TRUE)

#se$new.ident <- gsub("D\\d","D", se$orig.ident) #only for ESCC_TDN dataset
#https://github.com/satijalab/seurat/issues/3207#issuecomment-651862274
#meta <- se@meta.data
#tmp<-as.data.frame(meta%>%group_by(new.ident)%>%summarize(count=n())%>%ungroup())
#cells.use <- rownames(meta[meta$new.ident %in% tmp[tmp$count>k_filter,]$new.ident,])
#se <- subset(se,cells=cells.use)
#split_se <- SplitObject(se, split.by = "new.ident")

split_se <- SplitObject(se, split.by = "orig.ident")

#Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_se, nfeatures = nFeature)

# Prepare the SCT list object for integration
split_se <- PrepSCTIntegration(object.list = split_se, anchor.features = integ_features)

reference_datasets <- which(names(split_se) %in% refs)

# Find best buddies - can take a while to run
if(mode=="cca"){
        integ_anchors <- FindIntegrationAnchors(object.list = split_se, normalization.method = "SCT",
                anchor.features = integ_features, reference = reference_datasets, verbose = TRUE, k.filter=k_filter)
}else{
        split_se <- future_lapply(X = split_se, FUN = RunPCA, verbose=FALSE, features=integ_features)
        integ_anchors <- FindIntegrationAnchors(object.list = split_se, normalization.method = "SCT",
                anchor.features = integ_features, reference = reference_datasets, reduction = "rpca", verbose = TRUE, k.filter=k_filter)
}

# Integrate across samples and get a batch-corrected assay of all features (genes)
# https://github.com/satijalab/seurat/issues/2590#issuecomment-622463806
all_features <- Reduce(intersect, lapply(split_se, rownames))
se <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", features.to.integrate = all_features, verbose=T)
saveRDS(se, file=outfile)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
