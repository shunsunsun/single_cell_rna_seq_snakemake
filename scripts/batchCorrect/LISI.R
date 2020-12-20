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
#infile=snakemake@input[[1]] #GEJ_QCed_*_walktrap_clustered.rds
sce <- readRDS(file=infile)
meta <- as.data.frame(colData(sce))
meta <- meta[, colnames(meta)=='orig.ident' | grepl("snn", colnames(meta)) | grepl("_res.", colnames(meta))]
walktrap_res <- which(grepl("snn",colnames(meta)) & !grepl("_res",colnames(meta)))
colnames(meta)[walktrap_res]=gsub("snn","walktrap_snn",colnames(meta)[walktrap_res])
meta$orig.ident=sapply(strsplit(meta$orig.ident,"-"),"[",1)

for(r in reducedDimNames(sce)){
    outfile <- gsub("walktrap_clustered",paste0(r, "_lisi"),infile)
    if(!file.exists(outfile)){
        if(length(reducedDimNames(sce))>1){
                df=meta[,colnames(meta)=="orig.ident" | grepl(tolower(r),colnames(meta))]
        }
        rd <- reducedDim(sce, type=r)
        lisi <- compute_lisi(X=rd, meta_data=df, label_colnames='orig.ident', perplexity=40)
	res <- data.frame(BasedOn='all_cell',median=median(lisi[,1]),min=min(lisi[,1]),max=max(lisi[,1]),
		median_norm=(median(lisi[,1])-min(lisi[,1]))/(max(lisi[,1])-min(lisi[,1])))
	for(i in which(grepl("snn", colnames(df)) | grepl("_res.", colnames(df)))){
                clust_res <- df[,i]
                names(clust_res) = rownames(df)
                clust_lisi_median <- NULL
		clust_lisi_min <- NULL
		clust_lisi_max <- NULL
                clust_lisi_median_norm <- NULL
                for(c in unique(clust_res)){
                    clust_lisi <- lisi[names(clust_res[clust_res==c]),1]
		    md <- median(clust_lisi)
 		    mi <- min(clust_lisi)
	 	    mx <- max(clust_lisi)
                    clust_lisi_median <- c(clust_lisi_median, md)
		    clust_lisi_min <- c(clust_lisi_min, mi)
		    clust_lisi_max <- c(clust_lisi_max, mx)
                    clust_lisi_median_norm <- c(clust_lisi_median_norm, (md-mi)/(mx-mi))
                }
                res <- rbind(res, data.frame(BasedOn=colnames(df)[i],median=median(clust_lisi_median),min=median(clust_lisi_min),
			max=median(clust_lisi_max),median_norm=median(clust_lisi_median_norm)))
	}
        saveRDS(res, file=outfile)
    }else{
        print(paste0("Already done: ", basename(outfile)))
    }
}
