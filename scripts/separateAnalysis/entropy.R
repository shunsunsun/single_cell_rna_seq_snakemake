##Shannon entropy as proposed in https://github.com/cellgeni/batchbench/blob/master/bin/R_entropy.R
##A high batch entropy is preferrable

suppressPackageStartupMessages({
    library(future.apply)
    library(scran)
    library(RANN)
})
args  <- commandArgs(trailingOnly=T)
infile=args[1]
#infile=snakemake@input[[1]] #GEJ_QCed_*_walktrap_clustered.rds
useImmun=args[2]

sce <- readRDS(file=infile)
if(useImmun=='T'){
   meta <- colData(sce)
   sce <- sce[,grepl("-I",rownames(meta))]
}
meta <- as.data.frame(colData(sce))
meta <- meta[, colnames(meta)=='orig.ident' | grepl("snn", colnames(meta)) | grepl("_res.", colnames(meta))]
walktrap_res <- which(grepl("snn",colnames(meta)) & !grepl("_res",colnames(meta)))
colnames(meta)[walktrap_res]=gsub("snn","walktrap_snn",colnames(meta)[walktrap_res])
meta$orig.ident=sapply(strsplit(meta$orig.ident,"-"),"[",1)
nbatches <- length(unique(meta$orig.ident))

nworker=min(as.numeric(args[3]),length(availableWorkers()))
print(paste0("Use ",nworker, " workers"))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)

# entropy function
shannon_entropy <- function(x, batch_vector, N_batches) {
  x = x[-1]
  freq_batch = table(batch_vector[x])/length(batch_vector[x])
  freq_batch_positive = freq_batch[freq_batch > 0]
  return(-sum(freq_batch_positive * log(freq_batch_positive))/log(N_batches))
}


for(r in reducedDimNames(sce)){
    outfile <- gsub("walktrap_clustered",paste0(r, "_knn50Entropy"),infile)
    if(useImmun=='T'){
        outfile <- gsub("Entropy","EntropyImmun",outfile)
    }
    if(!file.exists(outfile)){
        if(length(reducedDimNames(sce))>1){
                df=meta[,colnames(meta)=="orig.ident" | grepl(tolower(r),colnames(meta))]
        }else{
		df=meta
	}
        rd <- reducedDim(sce, type=r)
	knn <- RANN::nn2(rd,k = 51)$nn.idx
	batch_entropy <- future_apply(knn, 1, future.seed=1129, FUN = function(x) {shannon_entropy (x, df$orig.ident, nbatches)})
	#names(batch_entropy)=rownames(df)
	res <- data.frame(BasedOn='all_cell',median=median(batch_entropy),min=min(batch_entropy),max=max(batch_entropy))
	for(i in which(grepl("snn", colnames(df)) | grepl("_res.", colnames(df)))){
                clust_res <- df[,i]
                #names(clust_res) = rownames(df)
                clust_entropy_median <- NULL
		clust_entropy_min <- NULL
		clust_entropy_max <- NULL
                for(c in unique(clust_res)){
                    clust_entropy <- batch_entropy[which(clust_res==c)]
                    clust_entropy_median <- c(clust_entropy_median, median(clust_entropy))
		    clust_entropy_min <- c(clust_entropy_min, min(clust_entropy))
		    clust_entropy_max <- c(clust_entropy_max, max(clust_entropy))
                }
                res <- rbind(res, data.frame(BasedOn=colnames(df)[i],median=median(clust_entropy_median),min=median(clust_entropy_min),
			max=median(clust_entropy_max)))
	}
        saveRDS(res, file=outfile)
    }else{
        print(paste0("Already done: ", basename(outfile)))
    }
}

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
