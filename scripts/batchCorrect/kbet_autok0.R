##kBET: https://github.com/theislab/kBET
##smaller kBET is preferrable

suppressPackageStartupMessages({
    library(future.apply)
    library(scran)
    library(kBET)
})
args  <- commandArgs(trailingOnly=T)
infile=args[1]
#infile=snakemake@input[[1]] #GEJ_QCed_*_walktrap_clustered.rds
outfile <- gsub("walktrap_clustered","kBET",infile)
logfile <- paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/log/",gsub("walktrap_clustered.rds","kbet.log",basename(infile)))
nworker <- as.numeric(args[2])
assayLabel <- args[3] #logcounts
set.seed(1129)

if(!file.exists(outfile)){
	sce <- readRDS(file=infile)
	df <- as.data.frame(colData(sce))
	df <- df[, colnames(df)=='orig.ident' | grepl("snn", colnames(df)) | grepl("_res.", colnames(df))]
	walktrap_res <- which(grepl("snn",colnames(df)) & !grepl("_res",colnames(df)))
	colnames(df)[walktrap_res]=gsub("snn","walktrap_snn",colnames(df)[walktrap_res])

	nworker=min(nworker,length(availableWorkers()))
	print(paste0("Use ",nworker, " workers"))
	plan("multiprocess", workers = nworker)
	options(future.globals.maxSize = 60*1024^3)
	
	kBET_res <- NULL
	for(i in which(grepl("snn", colnames(df)) | grepl("_res.", colnames(df)))){
		clustMethod <- colnames(df)[i]
		print(paste0("Evaluating ", clustMethod))
		clust_res <- df[,i]
		names(clust_res) <- rownames(df)
		kBET_tmp <- do.call("rbind", future_lapply(unique(clust_res), future.seed=TRUE, function(c,all,lab){
			sce_sub <- sce[,names(all)[all==c]]
			keep_feature <- rowSums(assay(sce_sub,lab) > 0) > 0
			sce_sub <- sce_sub[keep_feature,]
			data_sub <- t(assay(sce_sub,lab))
			batch_sub <- sce_sub$orig.ident
                        tryCatch({
				kBET_sub <- kBET(df=data_sub, batch=batch_sub, plot=FALSE)
				cat("kBET calculation finished",file=logfile,append=T,sep="\n")
	                        return(data.frame(cluster=c, kBET_obs=kBET_sub$summary$kBET.observed[1]))
			}, error=function(e){
				cat(sprintf("kBET calculation failed due to %s", conditionMessage(e)),file=logfile,append=T,sep="\n")
			})
		}, all=clust_res,lab=assayLabel))
		#save(kBET_tmp, file=paste0("temp_",clustMethod,".rda"))
		#kBET_tmp <- do.call("rbind",kBET_tmp)
		kBET_res <- rbind(kBET_res, data.frame(clustMethod=clustMethod, median_kBET=median(kBET_tmp$kBET_obs),
			min_kBET=min(kBET_tmp$kBET_obs), max_kBET=max(kBET_tmp$kBET_obs)))
	}
        saveRDS(kBET_res, file=outfile)
	options(future.globals.maxSize = 500*1024^2) #500M
	plan(sequential)
}else{
        print(paste0("Already done: ", basename(outfile)))
}
