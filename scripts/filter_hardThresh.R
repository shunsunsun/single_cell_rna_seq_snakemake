input <- snakemake@input[[1]]
filename <- unname(unlist(strsplit(basename(input),"_")))
cancer <- filename[1]
celltype <- filename[2]
outfile=snakemake@output[[1]]
outdir=dirname(outfile)
if(!dir.exists(outdir)){
        dir.create(outdir,recursive=T)
}

library(Seurat)
options(stringsAsFactors=F)

## Read input
seurat_obj <- readRDS(input)

####################
## Cell-level filtering

seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
seurat_obj$mitoCountRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")/100
# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = seurat_obj, subset= (nCount_RNA >= as.numeric(snakemake@params[["nUMI_lo"]])) & 
		(nCount_RNA >= as.numeric(snakemake@params[["nUMI_up"]])) &
		(nFeature_RNA >= as.numeric(snakemake@params[["nGene_lo"]])) & 
		(nFeature_RNA >= as.numeric(snakemake@params[["nGene_up"]])) &
		(log10GenesPerUMI > as.numeric(snakemake@params[["log10GenesPerUMI"]])) & 
		(mitoCountRatio < as.numeric(snakemake@params[["mitoCountRatio"]])))


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


saveRDS(filtered_seurat,file=outfile)
