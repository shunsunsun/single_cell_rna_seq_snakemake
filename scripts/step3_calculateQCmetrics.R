input <- snakemake@input[[1]]
output <- snakemake@output[[1]]
housekeeping <- snakemake@params[[1]]
cpu <- snakemake@threads

suppressPackageStartupMessages(library(Seurat))

## Read input
seurat_obj <- readRDS(input)

## Only need to do once: change orig.ident to sample ID
#seurat_obj$orig.ident <- sapply(rownames(seurat_obj@meta.data),function(x){
#	unlist(strsplit(x,"_"))[1]
#})
#saveRDS(seurat_obj,file=input)

########Add metadata: Seurat###################

seurat_obj$type <- ifelse(grepl("N",seurat_obj$orig.ident),"N","T")

## Add number of genes per UMI for each cell to metadata
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

## Compute the mitochodrial transcript percentage for each cell
## The PercentageFeatureSet() will take a pattern and search the gene identifiers.
## For each column (cell) it will take the sum of the counts slot for features belonging to the set, divide by the column sum for all features and multiply by 100.
seurat_obj$mitoCountRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")/100

## Compute the ribosomal transcript percentage for each cell
seurat_obj$riboCountRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^RP[SL][[:digit:]]")/100

## Compute the house-keeping transcript percentage for each cell
load(housekeeping)
hkgenes = as.character(Housekeeping_Genes$Gene.name)
hkgenes <- rownames(seurat_obj)[rownames(seurat_obj) %in% hkgenes]
hkgenes <- hkgenes[!is.na(hkgenes)]
seurat_obj$hbgeneCountRatio <-PercentageFeatureSet(object=seurat_obj, features=hkgenes)/100
###the following code achieves the same as the previous line
#ct <- GetAssayData(object = seurat_obj, slot = "counts")
#hkgeneRatio <- Matrix::colSums(ct[hkgenes,])/Matrix::colSums(ct)
#seurat_obj <- AddMetaData(seurat_obj, hkgeneRatio, col.name = "hkgeneCountRatio")

##Compute the count ratio of red blood cells marker genes
hbgenes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
hbgenes <- rownames(seurat_obj)[rownames(seurat_obj) %in% hbgenes]  
hbgenes <- hbgenes[!is.na(hbgenes)] 
seurat_obj$hbgeneCountRatio <-PercentageFeatureSet(seurat_obj, features=hbgenes)/100


suppressPackageStartupMessages(library(scater))

## Detect doublets
suppressPackageStartupMessages(library(scDblFinder))
suppressPackageStartupMessages(library(BiocParallel))
sce <- as.SingleCellExperiment(seurat_obj)
sce <- scDblFinder(sce, samples="ident", BPPARAM=MulticoreParam(cpu))

## Save the seurat object with added meta data
saveRDS(as.Seurat(sce),file=output)
