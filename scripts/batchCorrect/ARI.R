# Adjusted rand index test for overlap between batch and cluster labelings.The adjustment corrects for chance grouping between cluster elements. 
# This goes between 0 (completely dissimilar clustering) to 1 (identical clustering). We use 1-ARI since small values are good in this context
# Ref: https://davetang.org/muse/2017/09/21/adjusted-rand-index/
suppressPackageStartupMessages({
    library(Seurat)
    library(fossil)
    library(scater)
})
args  <- commandArgs(trailingOnly=T)
infile=args[1]
#infile=snakemake@input[[1]] #GEJ_QCed_*_walktrap_clustered.rds
outfile <- gsub("walktrap_clustered","ari",infile)
useImmun=args[2]
if(useImmun=='T'){
	outfile <- gsub("ari","ariImmun",outfile)
}

if(!file.exists(outfile)){
    sce <- readRDS(file=infile)
    if(useImmun=='T'){
	meta <- colData(sce)
	sce <- sce[,grepl("-I",rownames(meta))]
    }
    df <- as.data.frame(colData(sce))
    df <- df[, colnames(df)=="orig.ident" | grepl("snn", colnames(df)) | grepl("_res.", colnames(df))]
    samples <- unique(df$orig.ident)
    df$sample <- plyr::mapvalues(df$orig.ident, from = samples, to = 1:length(samples))
    walktrap_res <- which(grepl("snn",colnames(df)) & !grepl("_res",colnames(df)))
    colnames(df)[walktrap_res]=gsub("snn","walktrap_snn",colnames(df)[walktrap_res])
    res <- NULL
    for(i in which(grepl("snn", colnames(df)) | grepl("_res.", colnames(df)))){
        ari=adj.rand.index(as.numeric(df$sample), as.numeric(df[,i]))
        res=rbind(res, data.frame(clustMethod=colnames(df)[i], one_minus_ARI=1-ari))
    }
    saveRDS(res, file=outfile)
}else{
    print(paste0("Already done: ", basename(outfile)))
}
