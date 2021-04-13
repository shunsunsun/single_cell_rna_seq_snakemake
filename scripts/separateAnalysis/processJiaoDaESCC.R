#check /lustre4/gaog_pkuhpc/users/liny/ESCC_shanghai_jiaoda/GSE145370/readme.txt
suppressPackageStartupMessages({
	library(Seurat)
})
args  <- commandArgs(trailingOnly=T)
ref_file <- "/lustre4/gaog_pkuhpc/users/liny/ESCC_shanghai_jiaoda/GSE145370/ESCC_GSE145370.RDS"
celltype <- args[1]
cell_file <- gsub(".RDS",paste0("_",celltype,".RDS"),ref_file)

if(!file.exists(cell_file)){ #need to extract and process specific cell types
	ref <- readRDS(file=ref_file)
	if(celltype=="TCells"){
		clust_ids <- c(0:3,12,14,18)
	}
	cells.use <- rownames(ref@meta.data[ref@meta.data$seurat_clusters %in% clust_ids, ])
	ref <- subset(ref, cells=cells.use)
	ref$seurat_clusters <- NULL
	ref$RNA_snn_res.0.8 <- NULL
	saveRDS(ref,file=cell_file)
}
