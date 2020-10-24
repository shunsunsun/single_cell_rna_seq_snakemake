suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))

input <- snakemake@input[[1]] ## seurat object as stored as *.rds

nPC <- as.numeric(snakemake@params[["nPC"]])
pca_plot=snakemake@output[["pca"]]
tsne_plot=snakemake@output[["tsne"]]
umap_plot=snakemake@output[["umap"]]
qc_plot=snakemake@output[["qc"]]
outdir <- dirname(qc_plot)
if(!dir.exists(outdir)){
        dir.create(outdir,recursive=T)
}

nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))

##### Read input
se <- readRDS(input)
se <- ifelse(is(se, "Seurat"), se, as.Seurat(se))

plan("multiprocess", workers = nworker)
set.seed(1129)

##### Feature selection
##note: Selection of variable features is performed on @data or @counts, not @scale.data and therefore unaffected by ScaleData
metric=snakemake@params[["metric"]]
nFeature=snakemake@params[["topN"]]
if(metric=="vst"){
	se <- FindVariableFeatures(se, selection.method = "vst", nfeatures = nFeature)
}else if(metric=="deviance"){
	se <- 
}

##### Data scaling
if(!grepl("IntegBySeurat3", input, fixed=TRUE)){
        se=ScaleData(object=se)
}else{
	se=se[["integrated"]]
}

##### Dimension reduction by PCA, TSNE and UMAP embedding
se <- RunPCA(se,features=VariableFeatures(se))
se <- RunUMAP(seurat_obj, dims = 1:nPC)
tmp <- tryCatch({
	return(RunTSNE(seurat_obj, dims = 1:nPC))
}, error=function(x){
	return(RunTSNE(seurat_obj, dims = 1:nPC, perplexity = 5))
})
seurat_obj=tmp

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:nPC)
seurat_obj <- FindClusters(seurat_obj)

#jpeg(outfile, width = 18, height = 5, units="in", res=300)
#pca=DimPlot(seurat_obj, reduction = "pca", group.by= "Phase", split.by = "Phase")
#print(pca)
#dev.off()

#jpeg(outfile2, width = 8, height = 6, units="in", res=300)
#pca=DimPlot(seurat_obj, reduction = "pca", group.by= "Phase")
#print(pca)
#dev.off()

# Visualization
## Feature plot with QC metrics
metrics=c("nCount_RNA", "nFeature_RNA","subsets_Mito_percent","log10GenesPerUMI","S.Score","G2M.Score")
p=FeaturePlot(seurat_obj, reduction = "umap", features = metrics, ncol=3, pt.size = 1.2, order = TRUE, min.cutoff = 'q10', label = TRUE)
ggsave(file=qc_plot, plot = p, width = 20, height = 12, device="jpeg",units="in",dpi=300)

## PCA plot grouped and colored by cell cycle phases
#plot1 <- DimPlot(seurat_obj, reduction = "pca", group.by="Phase")
## PCA plot grouped and colored by clusters
#plot2 <- DimPlot(seurat_obj, reduction = "pca", label=T) + NoLegend()
#plot3 <- ElbowPlot(seurat_obj, ndims=50, reduction="pca")
#p = plot1 + plot2 + plot3
#ggsave(file=pca_plot, plot = p, width = 18, height = 5, device="jpeg",units="in",dpi=300)
## PCA plot split and colored by cell cycle phases
#plot1s <- DimPlot(seurat_obj, reduction = "pca", group.by="Phase", split.by="Phase")
#ggsave(file=gsub("pca","pca.split",pca_plot), plot = plot1s, width = 18, height = 5, device="jpeg",units="in",dpi=300)

#TSNE plot grouped and colored by clusters
#plot1 <- DimPlot(seurat_obj, reduction = "tsne", label=T) + NoLegend()
#TSNE plot grouped and colored by cell cycle phases
#plot2 <- DimPlot(seurat_obj, reduction = "tsne", group.by='Phase')
#g <- plot2 + plot1
#ggsave(file=tsne_plot, plot = g, width = 10, height = 5, device="jpeg",units="in",dpi=300)
##TSNE plot split and colored by cell cycle phases
#plot2s <- DimPlot(seurat_obj, reduction = "tsne", group.by='Phase',split.by="Phase")
#ggsave(file=gsub("tsne","tsne.split",tsne_plot), plot = plot2s, width = 18, height = 5, device="jpeg",units="in",dpi=300)

#UMAP grouped and colored by clusters
#plot1 <- DimPlot(seurat_obj, reduction = "umap", label=T) + NoLegend()
#UMAP grouped and colored by cell cycle phases
#plot2 <- DimPlot(seurat_obj, reduction = "umap", group.by='Phase')
#g <- plot2+plot1
#ggsave(file=umap_plot, plot = g, width = 10, height = 5, device="jpeg",units="in",dpi=300)
##UMAP split and colored by cell cycle phases
#plot2s <- DimPlot(seurat_obj, reduction = "umap", group.by='Phase',split.by="Phase")
#ggsave(file=gsub("umap","umap.split",umap_plot), plot = plot2s, width = 18, height = 5, device="jpeg",units="in",dpi=300)


