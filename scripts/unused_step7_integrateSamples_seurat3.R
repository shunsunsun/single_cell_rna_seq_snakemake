options(stringsAsFactors = F)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))

indir <- snakemake@params[["indir"]]
label <- snakemake@params[["infile"]]
outfile <- snakemake@output[["rds"]]
nPC <- as.numeric(snakemake@params[["nPC"]])
nFeature <- as.numeric(snakemake@params[["nFeature"]])

nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
#cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 25*1024^3) 

samples=scan(file=snakemake@input[["list"]],what="character")
normed_seurat <- future_lapply(X = samples, FUN = function(s){
        readRDS(file=paste0(indir,"/",s,".",label,".rds"))
})

#Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = normed_seurat, nfeatures = nFeature)

# Prepare the SCT list object for integration
normed_seurat <- PrepSCTIntegration(object.list = normed_seurat, anchor.features = integ_features)

# prevent small datasets from not having enough neighbors (k) to use when filtering anchors
# https://github.com/satijalab/seurat/issues/997#issuecomment-505207415
k_filter <- min(200, min(sapply(normed_seurat, ncol)))

# Find best buddies - can take a while to run
if(snakemake@params[["anchor"]]=="cca"){
	integ_anchors <- FindIntegrationAnchors(object.list = normed_seurat, normalization.method = "SCT", 
	anchor.features = integ_features, verbose = TRUE, dims=1:nPC, k.filter=k_filter)
}else{
	normed_seurat <- future_lapply(X = normed_seurat, FUN = RunPCA, verbose=FALSE, features=integ_features)
	integ_anchors <- FindIntegrationAnchors(object.list = normed_seurat, normalization.method = "SCT", 
	anchor.features = integ_features, reduction = "rpca", verbose = TRUE, dims=1:nPC, k.filter=k_filter)
}

# Integrate across samples and get a batch-corrected assay of all features (genes)
# https://github.com/satijalab/seurat/issues/2590#issuecomment-622463806
all_features <- Reduce(intersect, lapply(normed_seurat, rownames))
# same as: all_features <- lapply(normed_seurat, row.names) %>% Reduce(intersect,.)
normed_seurat <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", dims=1:nPC, features.to.integrate = all_features, verbose=T)
saveRDS(normed_seurat, file=outfile)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
