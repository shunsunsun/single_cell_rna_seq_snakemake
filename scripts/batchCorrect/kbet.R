##kBET: https://github.com/theislab/kBET

suppressPackageStartupMessages({
    library(future.apply)
    library(scran)
    library(kBET)
    library(FNN)
})
args  <- commandArgs(trailingOnly=T)
infile=args[1]
#infile=snakemake@input[[1]] #GEJ_QCed_*_walktrap_clustered.rds
outfile <- gsub("walktrap_clustered.rds","kBET.rda",infile)
set.seed(1129)

if(!file.exists(outfile)){
	sce <- readRDS(file=infile)
	df <- as.data.frame(colData(sce))
	df <- df[, colnames(df)=='orig.ident' | grepl("^snn", colnames(df)) | grepl("_res.", colnames(df))]
	colnames(df) <- gsub("^snn","walktrap_snn",colnames(df))     
	df$orig.ident=sapply(strsplit(df$orig.ident,"-"),"[",1)

	nworker=min(as.numeric(args[2]),length(availableWorkers()))
	print(paste0("Use ",nworker, " workers"))
	plan("multiprocess", workers = nworker)
	options(future.globals.maxSize = 60*1024^3)

	for(i in which(grepl("snn", colnames(df)) | grepl("_res.", colnames(df)))){
		clust_res <- df[,c(1,i)]
		colnames(clust_res)=c("batch","cluster")
		kBET_res <- setNames(future.lapply(X=unique(clust_res), future.seed=TRUE, function(c){
			batch_tmp <- clust_res[clust_res$cluster == c,]$batch
			data_tmp <- data[clusters == cluster_level,]
                        kBET_tmp <- kBET(df=data_tmp, batch=batch_tmp, plot=FALSE)
                        return(kBET_tmp)
		}), unique(clust_res))



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
