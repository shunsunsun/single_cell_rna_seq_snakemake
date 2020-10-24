suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(glmpca))

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
nPC <- as.numeric(snakemake@params[[1]])
randSeed <- 1129L

se <- readRDS(file=infile)
if(grepl("pipeComp",infile)){
	dims <- c(nPC,10,20)
	reducs <- c('pca','glmpca10','glmpca20')
}else if(grepl("outly",infile)){
	dims <- nPC
	reducs <- 'pca'
}

nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 10*1024^3) #10G
set.seed(randSeed)

for (i in 1:length(reducs)){
	clustName <- paste0(reducs[i],'Clust')
	tsneName <- paste0(reducs[i],'_tsne')
	umapName <- paste0(reducs[i],'_umap')
	if(!clustName %in% names(se)){
		se <- FindNeighbors(se, reduction=reducs[i], dims=1:dims[i], verbose=FALSE, graph.name=clustName, k.param = 50)
		se <- FindClusters(se, graph.name=clustName, resolution = seq(0.6,2,by=0.2), random.seed=randSeed, verbose=FALSE)
	}
	if(!tsneName %in% names(se)){
		se <- RunTSNE(se,reduction=reducs[i],dims=1:dims[i],seed.use=randSeed,reduction.name=paste0(reducs[i],'_tsne'),reduction.key = paste0(reducs[i],"tSNE_"))
	}
	if(!umapName %in% names(se)){
		se <- RunUMAP(se,reduction=reducs[i],dims=1:dims[i],seed.use=randSeed,reduction.name=paste0(reducs[i],'_umap'),reduction.key = paste0(reducs[i],"UMAP_"))
	}
}
saveRDS(se, outfile)
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
