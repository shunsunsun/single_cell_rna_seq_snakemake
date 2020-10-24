input <- snakemake@input[[1]]
output <- snakemake@output[[1]]
cpu <- snakemake@threads

suppressPackageStartupMessages(library(Seurat))

## Read input
seurat_obj <- readRDS(input)

## Only need to do once: change orig.ident to sample ID
#seurat_obj$orig.ident <- sapply(rownames(seurat_obj@meta.data),function(x){
#	unlist(strsplit(x,"_"))[1]
#})
#saveRDS(seurat_obj,file=input)

## Compute the mitochodrial transcript percentage for each cell
## The PercentageFeatureSet() will take a pattern and search the gene identifiers.
## For each column (cell) it will take the sum of the counts slot for features belonging to the set, divide by the column sum for all features and multiply by 100.
seurat_obj$mitoCountRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")/100

suppressPackageStartupMessages(library(scater))

## Detect doublets
suppressPackageStartupMessages(library(scDblFinder))
suppressPackageStartupMessages(library(BiocParallel))
sce <- as.SingleCellExperiment(seurat_obj)
sce <- scDblFinder(sce, samples="ident", BPPARAM=MulticoreParam(cpu))

## Remove doublets
sce=sce[,sce$scDblFinder.class != "doublet"]

sce <- addPerCellQC(sce, percent_top=c(20,50,100,200))
sce$log10_total_counts <- log10(sce$sum)
sce$log10_total_features <- log10(sce$detected)
sce$featcount_ratio <- sce$log10_total_counts/sce$log10_total_features

#' Returns the difference to the expected ratio of counts and number of features
#' @param df a cell metadata data.frame
#' @param do.plot Logical; whether to plot the count/feature relationship.
#' @param linear Logical; whether to model the relationship with a linear model (default TRUE), rather than a loess.
#' @return A vector of differences.
getFeatCountDist <- function(df, do.plot=FALSE, linear=TRUE){
        df <- as.data.frame(df)
        if(linear){
                mod <- lm(df$log10_total_features ~ df$log10_total_counts)
        }else{
                mod <- loess(df$log10_total_features ~ df$log10_total_counts)
        }
        pred <- predict(mod, newdata=data.frame(log10_total_counts=df$log10_total_counts))
        df$diff <- df$log10_total_features - pred
        if(do.plot){
                library(ggplot2)
                ggplot(df, aes(x=total_counts, y=total_features, colour=diff)) + geom_point() + geom_smooth(method = "loess", col="black")
        }
        df$diff
}
sce$featcount_dist <- getFeatCountDist(as.data.frame(colData(sce)))

##Save the doublet-removed SingleCellExperiment object
saveRDS(sce, file=output)
