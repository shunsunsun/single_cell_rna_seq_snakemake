input <- snakemake@input[[1]]
filename <- unname(unlist(strsplit(basename(input),"_")))
cancer <- filename[1]
celltype <- filename[2]
label <- paste0(cancer,"_",celltype)
outdir <- snakemake@params[["outdir"]]
if(!dir.exists(outdir)){
	dir.create(outdir,recursive=T)
}
housekeeping <- snakemake@params[["house"]]

library(Seurat)
library(ggplot2)

## Read input
seurat_obj <- readRDS(input)

## Only need to do once: change orig.ident to sample ID
#seurat_obj$orig.ident <- sapply(rownames(seurat_obj@meta.data),function(x){
#	unlist(strsplit(x,"_"))[1]
#})
#saveRDS(seurat_obj,file=input)

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
ct <- GetAssayData(object = seurat_obj, slot = "counts")
hkgeneRatio <- Matrix::colSums(ct[hkgenes,])/Matrix::colSums(ct)
seurat_obj <- AddMetaData(seurat_obj, hkgeneRatio, col.name = "hkgeneCountRatio")


############
## Seurat visualization

# Visualize QC metrics as a violin plot
jpeg(sprintf("%s/%s_VlnPlot_ribo_hk_perc.jpg", outdir, label), width = 18, height = 5, units="in", res=300)
vln <- VlnPlot(object = seurat_obj, features = c("riboCountRatio", "hkgeneCountRatio"), pt.size=0, ncol = 2, 
	group.by="orig.ident") + theme(legend.position="none") + xlab("")
print(vln)
dev.off()

jpeg(sprintf("%s/%s_VlnPlot_mito_perc.jpg", outdir, label), width = 10, height = 5, units="in", res=300)
vln <- VlnPlot(object = seurat_obj, features = "mitoCountRatio", pt.size=0, group.by="orig.ident") +
	theme(legend.position="none",axis.text.x=element_text(angle=90)) + scale_y_log10() + xlab("")
print(vln)
dev.off()

jpeg(sprintf("%s/%s_VlnPlot_nCount_5Kmax.jpg", outdir, label), width = 10, height = 5, units="in", res=300)
vln <- VlnPlot(object = seurat_obj, features = "nCount_RNA", pt.size=0, group.by="orig.ident", y.max=50000) + 
	theme(legend.position="none",axis.text.x=element_text(angle=90)) + scale_y_log10() + xlab("")
print(vln)
dev.off()

jpeg(sprintf("%s/%s_VlnPlot_nFeature.jpg", outdir, label), width = 10, height = 5, units="in", res=300)
vln <- VlnPlot(object = seurat_obj, features = "nFeature_RNA", pt.size=0, group.by="orig.ident") + 
	theme(legend.position="none",axis.text.x=element_text(angle=90)) + xlab("")
print(vln)
dev.off()


jpeg(sprintf("%s/%s_VlnPlot_nGenePerUMI.jpg", outdir, label), width = 10, height = 5, units="in", res=300)
vln <- VlnPlot(object = seurat_obj, features = "log10GenesPerUMI", pt.size=0, group.by="orig.ident") + 
	theme(legend.position="none",axis.text.x=element_text(angle=90)) + xlab("")
print(vln)
dev.off()

## Bar plot indicates the total number of cells per sample

jpeg(sprintf("%s/%s_BarPlot_nCellperSample.jpg",outdir,label), width = 10, height = 5, units="in", res=300)
p <- ggplot(data=seurat_obj@meta.data, aes(x=orig.ident, fill=orig.ident)) + geom_bar() + geom_hline(yintercept=snakemake@params[["nCell"]], 
	color = "black") + scale_y_log10() + theme(legend.position="none", axis.text.x=element_text(angle=90)) + xlab("") + ylab("nCell")
print(p)
dev.off()


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

jpeg(sprintf("%s/%s_ScattPlot_nCount_mitoRatio.jpg", outdir, label), width = 8, height = 6, units="in", res=300)
scatter <- FeatureScatter(object = seurat_obj, feature1 = "nCount_RNA", feature2 = "mitoCountRatio", pt.size=0.1) + 
	theme(legend.position="none") + xlab("nUMI") + scale_x_log10()
print(scatter)
dev.off()

jpeg(sprintf("%s/%s_ScattPlot_nCount_riboRatio.jpg", outdir, label), width = 8, height = 6, units="in", res=300)
scatter <- FeatureScatter(object = seurat_obj, feature1 = "nCount_RNA", feature2 = "riboCountRatio", pt.size=0.1) + 
	theme(legend.position="none") + xlab("nUMI") + scale_x_log10()
print(scatter)
dev.off()

jpeg(sprintf("%s/%s_ScattPlot_nCount_nFeature.jpg", outdir, label), width = 8, height = 6, units="in", res=300)
scatter <- FeatureScatter(object = seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1) + 
	theme(legend.position="none") + xlab("nUMI") + ylab("nGene") + scale_x_log10() + scale_y_log10()
print(scatter)
dev.off()

jpeg(sprintf("%s/%s_ScattPlot_nCount_nFeature_perSample.jpg", outdir, label), width = 30, height = 10, units="in", res=300) 
p=ggplot(data=seurat_obj@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=mitoCountRatio)) + geom_point(size=0.1) + xlab("nUMI") + ylab("nGene") +
	scale_colour_gradient(low = "gray90", high = "black") + stat_smooth() + scale_x_log10() + scale_y_log10() + theme(legend.position="bottom") +  
	facet_wrap(~orig.ident,ncol=9,scales="free")
print(p)
dev.off()
