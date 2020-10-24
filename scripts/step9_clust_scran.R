suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(BiocParallel))

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
nPC <- as.numeric(snakemake@params[[1]])
randSeed <- 1129L

se <- readRDS(file=infile)
sce <- as.SingleCellExperiment(se) #the graph information in the seurat object is not retained
if(grepl("pipeComp",infile)){
	dims <- c(nPC,10,20)
	reducs <- c('pca','glmpca10','glmpca20')
}else if(grepl("outly",infile)){
	dims <- nPC
	reducs <- 'pca'
}

set.seed(randSeed)
graph <- list()
for (i in 1:length(reducs)){
	for(k in c(10,30,50)){
		#MulticoreParam() by default uses all cores available as determined by detectCores
		g <- buildSNNGraph(sce, k=k, use.dimred = toupper(reducs[i]),BPPARAM=MulticoreParam())
		#graph[[paste0(reducs[i],"Clust_k.",k)]] <- g
		#sce[[paste0(reducs[i],"Clust_k.",k)]] <- igraph::cluster_walktrap(g)$membership
		graph[[paste0(reducs[i],"_k.",k)]] <- g
		sce[[paste0(reducs[i],"_k.",k)]] <- igraph::cluster_walktrap(g)$membership
	}
}
save(sce, graph, file=outfile)
