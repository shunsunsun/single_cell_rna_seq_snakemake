##kBET: https://github.com/theislab/kBET
##smaller kBET is preferrable

suppressPackageStartupMessages({
    library(foreach)
    library(doParallel)
    library(scran)
    library(kBET)
})
args  <- commandArgs(trailingOnly=T)
infile=args[1]
outfile=gsub(".rds","_kBET.rds",infile)
sce <- readRDS(file=infile)
sce <- as.SingleCellExperiment(sce)
meta <- colData(sce)
nworker <- as.numeric(args[2])
assayLabel <- args[3] #logcounts
set.seed(1129)

if(!file.exists(outfile)){
	nworker=min(nworker,detectCores())
	print(paste0("Use ",nworker, " workers"))
	registerDoParallel(nworker)
	
	clust_res <- meta[,grepl("_res.", colnames(meta))]
	names(clust_res) <- rownames(meta)
	uniq_c <- unique(clust_res)
	kBET_res <- foreach(j=1:length(uniq_c), .combine=rbind) %dopar% {
		sce_sub <- sce[,names(clust_res)[clust_res==uniq_c[j]]]
		keep_feature <- rowSums(assay(sce_sub,assayLabel) > 0) > 0
		sce_sub <- sce_sub[keep_feature,]
		data_sub <- t(assay(sce_sub,assayLabel))
		clust_size <- nrow(data_sub)
		if(floor(25*clust_size/100)<10){
			return(NULL)
		}
		batch_sub <- sce_sub$orig.ident
		kBET_sub <- c()
		for(p in seq(from = 5, to = 25, by = 5)){
			k0_val <- floor(p*nrow(data_sub)/100)
			if(k0_val<10){
				print(paste0("Skip: k0=",k0_val," is too small"))
				next
			}
                       	tryCatch({
				kBET_sub_p <- kBET(df=data_sub, k0=k0_val, batch=batch_sub, plot=FALSE, heuristic=FALSE)
	                       	kBET_sub <- c(kBET_sub,kBET_sub_p$summary$kBET.observed[1])
				print(paste0("kBET calculation with k0=",k0_val," finished"))
			}, error=function(e){
				print(paste0("kBET calculation with k0=",k0_val," failed due to ", conditionMessage(e)))
			})
		}
		data.frame(cluster=uniq_c[j], kBET_obs=median(kBET_sub))
	}
        saveRDS(kBET_res, file=outfile)
	stopImplicitCluster()
}else{
        print(paste0("Already done: ", basename(outfile)))
}
