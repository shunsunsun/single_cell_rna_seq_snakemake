suppressPackageStartupMessages({
	library(tidyverse)
	library(future)
	library(future.apply)
})

args  <- commandArgs(trailingOnly=T)
#rdss<- unname(unlist(srsplit(args[1],',')))
rdss <- paste0("../../data/ESCC_QCed_sctNorm_BatchHmy_clustStab/subsample/subsample_k_100_PC_50_round_", 0:19, ".rds")
output <- "../../data/ESCC_QCed_sctNorm_BatchHmy_clustStab/gather_subsample.rds"
resolution <- c(0.2,0.4,0.6)
#resolution <- as.numeric(unname(unlist(strsplit(as.character(args[3]),','))))
nworker=min(as.numeric(args[1]),length(availableWorkers()))

plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
set.seed(1129L)

get_df<- function(rds, resolution=NULL){
	res <- readRDS(rds)
	#res$pc=as.character(ress$pc)
	#res$round=as.character(res$round)
	#res$k_param=as.character(res$k_param)
	#res$resolution=as.character(res$resolution)
	if(length(resolution) > 0){
	#	resCol <- colnames(res[[]])[grepl("_res.",colnames(res[[]]))]
	#	unTgtCol <- colnames(res[[]])[!grepl(paste0("_res.", paste(resolution,collapse="|")), colnames(res[[]]))]
	#	col2rm <- intersect(resCol, unTgtCol)
	#	for(c in col2rm){
	#		res[[c]] <- NULL
	#	}
		res <- res[res$resolution %in% resolution, ]
	}
	return(res)
}

dat.list<- future_lapply(rdss, FUN=get_df, resolution=resolution)
gather_idents<- do.call(bind_rows, dat.list)
saveRDS(gather_idents, file = output)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)

