##Ref: https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
##Ref: https://satijalab.org/seurat/articles/sctransform_vignette.html

options(stringsAsFactors = F)
suppressPackageStartupMessages({
        library(Seurat)
	library(dplyr)
        library(future)
        library(future.apply)
})

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
outfile <- args[2]
nFeature <- 3000

nworker <- min(as.numeric(args[3]),length(availableWorkers()))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
set.seed(1129L)

se <- readRDS(file=infile)
DefaultAssay(se) <- "RNA"
se <- DietSeurat(se, assays="RNA")

#For ESCC_TDN dataset
if(grepl("TDN",infile)){
	se$patient <- gsub("D\\d","D",se$orig.ident)
	se$patient <- gsub("[NTD]","",se$patient)
}else{	#For ESCC jiaoda
	se$orig.ident <- gsub("_filtered_feature_bc_matrix","",se$orig.ident)
}

split_seurat <- SplitObject(se, split.by = "orig.ident")	

for(i in 1:length(split_seurat)){
   split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
   split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]], selection.method = "vst", nfeatures = nFeature)	
}

#TBD

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
