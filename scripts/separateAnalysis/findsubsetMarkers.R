suppressPackageStartupMessages({
        library(Seurat)
        library(future)
        library(tidyverse)
})

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
sctransform <- args[2] #T or F
assay_used <- ifelse(sctransform=="T","SCT","RNA") #RNA is recommended; I only find SCT useful for integrated datasets
outfile <- args[3]

subsets <- unlist(strsplit(args[4],","))
if(length(subsets)<2){
	print("Error: must have at least two subsets for comparison")
}
nworker <- min(as.numeric(args[5]),length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3) #60G
#randSeed <- 1129L
#set.seed(randSeed)

se <- readRDS(file=infile)
clustRes <- colnames(se@meta.data)[grepl("_res\\.", colnames(se@meta.data))]
Idents(se) <- clustRes
res <- NULL
for(s in subsets){
	markers <- FindMarkers(object = se, assay=assay_used, slot="data", ident.1 = s, ident.2 = setdiff(subsets, s), min.pct = 0.25)
	df=cbind(cluster=s, markers)
	res <- rbind(res, df)
}
saveRDS(res, file=outfile)
