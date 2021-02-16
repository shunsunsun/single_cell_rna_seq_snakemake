suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(BiocParallel))

args  <- commandArgs(trailingOnly=T)

infile <- args[1]
outfile <- args[2]
reducs <- args[3]
k <- as.numeric(args[4]) #50, 100, 150
randSeed <- 1129L

se <- readRDS(file=infile)
sce <- as.SingleCellExperiment(se) #the graph information in the seurat object is not retained
set.seed(randSeed)
#MulticoreParam() by default uses all cores available as determined by detectCores
g <- buildSNNGraph(sce, k=k, use.dimred = toupper(reducs),BPPARAM=MulticoreParam())
sce[[paste0(reducs,"_k.",k)]] <- igraph::cluster_walktrap(g)$membership
saveRDS(sce, file=outfile)
