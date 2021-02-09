suppressPackageStartupMessages({
	library(tidyverse)
	library(Seurat)
	library(future.apply)
})

args  <- commandArgs(trailingOnly=T)
rdss<- args[1]
output <- args[2]
resolution <- as.numeric(unname(unlist(strsplit(as.character(args[3]),','))))
nworker=min(as.numeric(args[4]),length(availableWorkers()))

plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
set.seed(1129)

get_clusts<- function(rds, resolution=NULL){
        se <- readRDS(rds)
        k<- gsub("full_sample_k_([0-9]+)_PC_([0-9]+).rds", "\\1", basename(rds))
        p<- gsub("full_sample_k_([0-9]+)_PC_([0-9]+).rds", "\\2", basename(rds))	
	tmp <- NULL
	if(length(resolution)==0){
		tgtCol <- colnames(se[[]])[grepl("_res.",colnames(se[[]]))]
	}else{
		tgtCol <- colnames(se[[]])[grepl(paste(sprintf("_res.%s",resolution),collapse="|"), colnames(se[[]]))]
	}
	for (c in tgtCol){
	    r <- unname(unlist(strsplit(c,'res\\.')))[2]
            tmp <- rbind(tmp, tibble::tibble(pc = p, resolution = r, k_param = k, original_ident_full = list(se[[c]])))		
        }
	return(tmp)
}

dat.list <- future_lapply(rdss, FUN=get_clusts, resolution=resolution)
gather_clusts<- do.call(bind_rows, dat.list)
saveRDS(gather_clusts, file = snakemake@output[[1]])

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
