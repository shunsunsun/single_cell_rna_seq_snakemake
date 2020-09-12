input <- snakemake@input[[1]]
filename <- unname(unlist(strsplit(basename(input),"_")))
cancer <- filename[1]
celltype <- filename[2]
outdir=dirname(snakemake@output[["rds"]])
if(!dir.exists(outdir)){
        dir.create(outdir,recursive=T)
}

library(Seurat)
options(stringsAsFactors=F)

## Read input
seurat_obj <- readRDS(input)

####################
## Cell-level filtering

library(scater)
sce=as.SingleCellExperiment(seurat_obj)
is.mito <- grepl("^MT-",rownames(sce))
sce=addPerCellQC(sce,subsets=list(Mito=is.mito))
if(snakemake@params[["all"]]){
	qc.lib <- isOutlier(sce$sum, log=TRUE, type="lower", batch=sce$orig.ident)
	qc.nexpr <- isOutlier(sce$detected, log=TRUE, type="lower", batch=sce$orig.ident)
	qc.mito <- filtered_seuratisOutlier(sce$subsets_Mito_percent, type="higher", batch=sce$orig.ident)
} else if(snakemake@params[["lowqual_list"]] != "NA"){
	qc.lib <- isOutlier(sce$sum, log=TRUE, type="lower", batch=sce$orig.ident)
	qc.nexpr <- isOutlier(sce$detected, log=TRUE, type="lower", batch=sce$orig.ident)
	lowqual_samples=read.table(file=snakemake@params[["lowqual_list"]],header=T)
	lowqual_samples=lowqual_samples[lowqual_samples$Cancer==cancer & lowqual_samples$Cell_type==celltype,]$Sample
	if(!is.null(lowqual_samples) && length(lowqual_samples)>0){
		qc.mito <- isOutlier(sce$subsets_Mito_percent, type="higher", batch=sce$orig.ident, subset=!(sce$orig.ident %in% lowqual_samples))
	}else{
		qc.mito <- isOutlier(sce$subsets_Mito_percent, type="higher", batch=sce$orig.ident)
	}
}
sce=sce[,!(qc.lib | qc.nexpr | qc.mito)]
filtered_seurat <- as.Seurat(sce)

# write stats
stats <- data.frame(ByLibSize=sum(qc.lib), ByFeature=sum(qc.nexpr), ByMito=sum(qc.mito), Remaining=ncol(sce))
write.table(stats, file=snakemake@output[["stats"]], sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)


###############
## Gene-level filtering

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than X TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= snakemake@params[["nCell"]]

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

saveRDS(filtered_seurat,file=snakemake@output[["rds"]])
