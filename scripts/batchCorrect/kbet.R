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
#logfile <- paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/log/",gsub("walktrap_clustered.rds","kbet.log",basename(infile)))
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
		#kBET_tmp <- NULL
		#all <- clust_res
		#lab <- assayLabel
		#for(c in unique(clust_res)){
		kBET_tmp <- do.call("rbind", future_lapply(unique(clust_res), future.seed=TRUE, function(c,all,lab){
			sce_sub <- sce[,names(all)[all==c]]
			keep_feature <- rowSums(assay(sce_sub,lab) > 0) > 0
			sce_sub <- sce_sub[keep_feature,]
			data_sub <- t(assay(sce_sub,lab))
			batch_sub <- sce_sub$orig.ident
			kBET_sub <- c()
			for(p in seq(from = 5, to = 25, by = 5)){
				k0_val <- floor(p*nrow(data_sub)/100)
                        	tryCatch({
					kBET_sub_p <- kBET(df=data_sub, k0=k0_val, batch=batch_sub, plot=FALSE, heuristic=FALSE)
					print(paste0("kBET calculation with k0=",k0_val," finished"))
	                        	kBET_sub=c(kBET_sub,kBET_sub_p$summary$kBET.observed[1])
				}, error=function(e){
					print(paste0("kBET calculation with k0=",k0_val," failed due to ", conditionMessage(e)))
					kBET_sub=c(kBET_sub,1)
				})
			}
			#print(kBET_sub)
			#kBET_tmp <- rbind(kBET_tmp, data.frame(cluster=c, kBET_obs=median(kBET_sub)))	
			return(data.frame(cluster=c, kBET_obs=median(kBET_sub)))
		}, all=clust_res,lab=assayLabel))
		#save(kBET_tmp, file=gsub(".rds", paste0(clustMethod,".rda"), outfile))
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
