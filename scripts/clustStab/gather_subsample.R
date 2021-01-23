suppressPackageStartupMessages({
	library(tidyverse)
	library(future.apply)
})

#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

rdss<- snakemake@input[["rds"]]
resolution <- as.numeric(unname(unlist(strsplit(as.character(snakemake@params[["res"]]),','))))

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
saveRDS(gather_idents, file = snakemake@output[[1]])
