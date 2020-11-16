options(stringsAsFactors = F)
suppressPackageStartupMessages({
	library(sctransform)
	library(Seurat)
	library(future)
	library(future.apply)
})

infile <- snakemake@input[[1]]
outfile <- snakemake@output[["rds"]]
nFeature <- as.numeric(snakemake@params[["nFeature"]])

nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
#cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
#options(future.globals.maxSize = 50*1024^3) #work for the GEJ but not the ESCC dataset
options(future.globals.maxSize = 70*1024^3)
set.seed(1129)

se <- readRDS(file=infile)
# https://github.com/satijalab/seurat/issues/2814#issuecomment-612103616
# After SCTransform, do not run FindVariableFeatures, which is designed for the LogNormal data.
# VariableFeatures(se[["SCT"]]) is rownames(se[["SCT"]]@scale.data)
# You can change the number of variables you need using the parameter, variable.features.n in the SCTransform

se@meta.data$orig.ident=sapply(strsplit(se@meta.data$orig.ident,"-"),"[",1)
split_se <- SplitObject(se, split.by = "orig.ident")

#Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_se, nfeatures = nFeature)

# Prepare the SCT list object for integration
split_se <- PrepSCTIntegration(object.list = split_se, anchor.features = integ_features)

if(grepl("ESCC", infile)){
	refs <- c("P31","P79")
}
if(grepl("GEJ", infile)){
	refs <- c("P29","P81")
}
reference_datasets <- which(names(split_se) %in% refs)

# prevent small datasets from not having enough neighbors (k) to use when filtering anchors
# https://github.com/satijalab/seurat/issues/997#issuecomment-505207415
k_filter <- min(200, min(sapply(split_se, ncol)))

# Find best buddies - can take a while to run
if(snakemake@params[["anchor"]]=="cca"){
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
