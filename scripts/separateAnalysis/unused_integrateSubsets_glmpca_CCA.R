suppressPackageStartupMessages({
        library(Seurat)
	library(SeuratWrappers)
	library(glmpca)
        library(future)
        library(future.apply)
})

args  <- commandArgs(trailingOnly=T)
obj1_file <- args[1]
obj1_reduc <- args[2] #e.g., glmpca10
obj2_file <- args[3]
obj2_reduc <- args[4]
resfile <- args[5]
nFeature <- as.numeric(args[6]) #2000
reduc <- args[7] ##cca or rpca
nworker <- min(as.numeric(args[8]),length(availableWorkers()))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
set.seed(1129L)

preprocess <- function(file){
	se <- readRDS(file=file)
    	regVars=c('mitoCountRatio', 'nFeature_RNA', 'nCount_RNA')
	se <- SCTransform(se, variable.features.n=3000, vars.to.regress = c(regVars, 'cc_difference'))		
	m <- GetAssayData(se, slot = "counts", assay = "RNA")
	devs <- scry::devianceFeatureSelection(m)
	VariableFeatures(se, assay="RNA") <- rownames(se)[order(devs, decreasing = TRUE)]
	return(se)
}

obj.list <- list()
obj.list[[1]] <- preprocess(obj1_file)
obj.list[[2]] <- preprocess(obj2_file)

#Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = nFeature)

# Prepare the list object for integration
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = obj.list, dims=1:30, normalization.method = "SCT",
                anchor.features = integ_features, reduction = reduc, verbose = TRUE)

all_features <- Reduce(intersect, lapply(obj.list, rownames))

#se <- IntegrateData(anchorset = integ_anchors, normalization.method="SCT", weight.reduction=c(obj1_reduc,obj2_reduc), features.to.integrate = all_features, verbose=TRUE)

if(!dir.exists(dirname(resfile))){
        dir.create(dirname(resfile),recursive=T)
}
save(all_features, integ_anchors, obj1_reduc, obj2_reduc, file=gsub("rds","rda",resfile))
#saveRDS(se, file=resfile)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
