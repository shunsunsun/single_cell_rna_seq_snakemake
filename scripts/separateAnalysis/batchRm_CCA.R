options(stringsAsFactors = F)
suppressPackageStartupMessages({
	library(sctransform)
	library(Seurat)
	library(dplyr)
	library(future)
	library(future.apply)
})

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
outfile <- args[2]
nFeature <- 2000
min_ncell <- as.numeric(args[3])
k.filter <- as.numeric(args[4])
method <- args[5] # cca or rpca
nworker <- min(as.numeric(args[6]),length(availableWorkers()))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 70*1024^3)
set.seed(1129)

se <- readRDS(file=infile)
#se <- DietSeurat(se, assays=c("RNA","SCT"), dimreducs="pca")
se <- DietSeurat(se, assays=c("RNA","SCT"))

se@meta.data$orig.ident=sapply(strsplit(se@meta.data$orig.ident,"-"),"[",1)
meta <- se@meta.data
tmp<-as.data.frame(meta%>%group_by(orig.ident)%>%summarize(count=n())%>%ungroup())
cells.use <- rownames(meta[meta$orig.ident %in% tmp[tmp$count>k.filter,]$orig.ident,])
se <- subset(se,cells=cells.use)

split_se <- SplitObject(se, split.by = "orig.ident")
reference_datasets <- which(names(split_se) %in% tmp[tmp$count>min_ncell,]$orig.ident)

#Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_se, nfeatures = nFeature)

# Prepare the SCT list object for integration
split_se <- PrepSCTIntegration(object.list = split_se, anchor.features = integ_features)

# Find best buddies - can take a while to run
if(method=="cca"){
	integ_anchors <- FindIntegrationAnchors(object.list = split_se, normalization.method = "SCT", 
		k.filter = k.filter, anchor.features = integ_features, reference = reference_datasets, verbose = TRUE)
}else{
	split_se <- future_lapply(X = split_se, FUN = RunPCA, verbose=FALSE, features=integ_features)
	integ_anchors <- FindIntegrationAnchors(object.list = split_se, normalization.method = "SCT", k.filter=k.filter,
		anchor.features = integ_features, reference = reference_datasets, reduction = "rpca", verbose = TRUE)
}

# Integrate across samples and get a batch-corrected assay of all features (genes)
# https://github.com/satijalab/seurat/issues/2590#issuecomment-622463806
all_features <- Reduce(intersect, lapply(split_se, rownames))
se <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", features.to.integrate = all_features, verbose=T)
saveRDS(se, file=outfile)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
