suppressPackageStartupMessages({
	library(tidyverse)
	library(Seurat)
)}

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


rdss<- snakemake@input[["rds"]]

get_clusts<- function(rds){
        x<- readRDS(rds)
        k<- gsub("full_sample_k_([0-9]+)_PC_([0-9]+).rds", "\\1", basename(rds))
        pc.use<- gsub("full_sample_k_([0-9]+)_PC_([0-9]+).rds", "\\2", basename(rds))	
	tmp_list <- list()
	for (c in colnames(se[[]])[grepl("_res.",colnames(se[[]])]){
	    r <- unname(unlist(strsplit(c,'.')[2]
            tmp_list.append(tibble::tibble(pc = p, resolution = r, k_param = k, original_ident_full = list(se[[c]])))
        }
        return(do.call(bind_rows, tmp_list))
}

dat.list <- lapply(rdss, get_clusts)
gather_clusts<- do.call(bind_rows, dat.list)
saveRDS(gather_clusts, file = snakemake@output[[1]])
