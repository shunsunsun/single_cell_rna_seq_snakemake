suppressPackageStartupMessages({
	library(Seurat)
	#library(BiocParallel)
	library(ggplot2)
	library(patchwork)
})

args  <- commandArgs(trailingOnly=T)
n_cpu <- as.numeric(args[1])
infile <- args[2] #absolute path
outPrefix <- gsub(".rds","_eval",infile)

se <- readRDS(file=infile)
df <- se@meta.data
df$orig.ident=sapply(strsplit(df$orig.ident,"-"),"[",1)

# Adjusted rand index test for overlap between batches (samples) and cluster labelings
# This goes between 0 (completely dissimilar clustering) to 1 (identical clustering)
# The adjustment corrects for chance grouping between cluster elements
# Ref: https://davetang.org/muse/2017/09/21/adjusted-rand-index/
res <- NULL
for(i in which(grepl("res.",colnames(df)))){
	ari <- dplyr::select(df, orig.ident, colnames(df)[i])
	colnames(ari)=c("batch","cluster")
	samples <- unique(ari$batch)
	ari$batch <- plyr::mapvalues(ari$batch, from = samples, to = 1:legnth(samples))
	score <- adj.rand.index(as.numeric(ari$batch), as.numeric(ari$cluster))
	res <- rbind(res, data.frame(ari=score, clustering=colnames(df)[i]))
}
saveRDS(res, file=paste0(outPrefix,"_ARI.rds"))
