options(stringsAsFactors = F)
suppressPackageStartupMessages(library(Seurat))

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
col2rm <- args[2] #e.g., SCT_snn_res.0.05

se <- readRDS(file=infile)
se$seurat_clusters <- NULL
se[[col2rm]] <- NULL
saveRDS(se, file=infile)
