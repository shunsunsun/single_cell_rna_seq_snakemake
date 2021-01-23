#local inverse Simpson’s index (LISI) measure the local batch distribution, based on local neighbors chosen on a preselected perplexity
#perplexity is the effective number of neighbors around each cell. Since we have 41 samples/batches, I use perplexity=40
#We use LISI to see if cell from different batches are well mixed globally and in each cluster

suppressPackageStartupMessages({
    library(Seurat)
    library(scater)
    library(lisi)
})
args  <- commandArgs(trailingOnly=T)
infile=args[1]
outfile=gsub(".rds","_lisi.rds",infile)
sce <- readRDS(file=infile)
sce <- as.SingleCellExperiment(sce)
meta <- colData(sce)

if(!file.exists(outfile)){
        rd <- reducedDim(sce, type="PCA")
        lisi <- compute_lisi(X=rd, meta_data=meta, label_colnames='orig.ident', perplexity=40)
	res <- data.frame(BasedOn='all_cell',median=median(lisi[,1]),min=min(lisi[,1]),max=max(lisi[,1]),
		median_norm=(median(lisi[,1])-min(lisi[,1]))/(max(lisi[,1])-min(lisi[,1])))
        clust_res <- meta[,grepl("_res.", colnames(meta))]
        names(clust_res) = rownames(meta)
        for(c in unique(clust_res)){
            clust_lisi <- lisi[names(clust_res[clust_res==c]),1]
	    md <- median(clust_lisi)
 	    mi <- min(clust_lisi)
	    mx <- max(clust_lisi)
            res <- rbind(res, data.frame(BasedOn=paste0("Cluster",c),median=md,min=mi,
		max=mx,median_norm=(md-mi)/(mx-mi)))
	}
        saveRDS(res, file=outfile)
}else{
        print(paste0("Already done: ", basename(outfile)))
}
