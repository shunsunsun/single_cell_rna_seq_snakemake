input <- snakemake@input[[1]]
filename <- unname(unlist(strsplit(basename(input),"_")))
cancer <- filename[1]
celltype <- filename[2]
label <- paste0(cancer,"_",celltype)
outdir <- snakemake@params[["outdir"]]
if(!dir.exists(outdir)){
	dir.create(outdir,recursive=T)
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))

## Read input
seurat_obj <- readRDS(input)

############
## Seurat visualization

# Visualize QC metrics as a violin plot

if(snakemake@params[["vln_ribo"]]){
	jpeg(sprintf("%s/%s_VlnPlot_ribo.jpg", outdir, label), width = 18, height = 5, units="in", res=300)
	vln <- VlnPlot(object = seurat_obj, features = "riboCountRatio", pt.size=0, ncol = 2, 
		group.by="orig.ident") + theme(legend.position="none") + xlab("")
	print(vln)
	dev.off()
}

if(snakemake@params[["vln_hk"]]){
	jpeg(sprintf("%s/%s_VlnPlot_hk.jpg", outdir, label), width = 18, height = 5, units="in", res=300)
	vln <- VlnPlot(object = seurat_obj, features = "hkgeneCountRatio", pt.size=0, ncol = 2,
		group.by="orig.ident") + theme(legend.position="none") + xlab("")
	print(vln)
	dev.off()
}

if(snakemake@params[["vln_mito"]]){
	jpeg(sprintf("%s/%s_VlnPlot_mito.jpg", outdir, label), width = 10, height = 5, units="in", res=300)
	vln <- VlnPlot(object = seurat_obj, features = "mitoCountRatio", pt.size=0, group.by="orig.ident") +
		theme(legend.position="none",axis.text.x=element_text(angle=90)) + scale_y_log10() + xlab("")
	print(vln)
	dev.off()
}

if(snakemake@params[["vln_hb"]]){
	jpeg(sprintf("%s/%s_VlnPlot_hb.jpg", outdir, label), width = 10, height = 5, units="in", res=300)
	vln <- VlnPlot(object = seurat_obj, features = "hbgeneCountRatio", pt.size=0, group.by="orig.ident") +
		theme(legend.position="none",axis.text.x=element_text(angle=90)) + scale_y_log10() + xlab("")
	print(vln)
	dev.off()
}

if(snakemake@params[["vln_nCount"]]){
	jpeg(sprintf("%s/%s_VlnPlot_nCount_5Kmax.jpg", outdir, label), width = 10, height = 5, units="in", res=300)
	vln <- VlnPlot(object = seurat_obj, features = "nCount_RNA", pt.size=0, group.by="orig.ident", y.max=50000) + 
		theme(legend.position="none",axis.text.x=element_text(angle=90)) + scale_y_log10() + xlab("")
	print(vln)
	dev.off()
}

if(snakemake@params[["vln_nFeature"]]){
	jpeg(sprintf("%s/%s_VlnPlot_nFeature.jpg", outdir, label), width = 10, height = 5, units="in", res=300)
	vln <- VlnPlot(object = seurat_obj, features = "nFeature_RNA", pt.size=0, group.by="orig.ident") + 
		theme(legend.position="none",axis.text.x=element_text(angle=90)) + xlab("")
	print(vln)
	dev.off()
}

if(snakemake@params[["vln_complexity"]]){
	jpeg(sprintf("%s/%s_VlnPlot_nGenePerUMI.jpg", outdir, label), width = 10, height = 5, units="in", res=300)
	vln <- VlnPlot(object = seurat_obj, features = "log10GenesPerUMI", pt.size=0, group.by="orig.ident") + 
		theme(legend.position="none",axis.text.x=element_text(angle=90)) + xlab("")
	print(vln)
	dev.off()
}

## Bar plot indicates the total number of cells per sample

if(snakemake@params[["bar_nCell"]]){
	jpeg(sprintf("%s/%s_BarPlot_nCellperSample.jpg",outdir,label), width = 10, height = 5, units="in", res=300)
	p <- ggplot(data=seurat_obj@meta.data, aes(x=orig.ident, fill=orig.ident)) + geom_bar() + geom_hline(yintercept=snakemake@params[["nCell"]], 
		color = "black") + scale_y_log10() + theme(legend.position="none", axis.text.x=element_text(angle=90)) + xlab("") + ylab("nCell")
	print(p)
	dev.off()
}

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

if(snakemake@params[["scatter_nCount_mito"]]){
	jpeg(sprintf("%s/%s_ScattPlot_nCount_mito.jpg", outdir, label), width = 8, height = 6, units="in", res=300)
	scatter <- FeatureScatter(object = seurat_obj, feature1 = "nCount_RNA", feature2 = "mitoCountRatio", pt.size=0.1) + 
		theme(legend.position="none") + xlab("nUMI") + scale_x_log10()
	print(scatter)
	dev.off()
}


#jpeg(sprintf("%s/%s_ScattPlot_nCount_ribo.jpg", outdir, label), width = 8, height = 6, units="in", res=300)
#scatter <- FeatureScatter(object = seurat_obj, feature1 = "nCount_RNA", feature2 = "riboCountRatio", pt.size=0.1) +
#	theme(legend.position="none") + xlab("nUMI") + scale_x_log10()
#print(scatter)
#dev.off()

if(snakemake@params[["scatter_nCount_nFeature_colSample"]]){
	jpeg(sprintf("%s/%s_ScattPlot_nCount_nFeature.jpg", outdir, label), width = 8, height = 6, units="in", res=300)
	scatter <- FeatureScatter(object = seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1) + 
		theme(legend.position="none") + xlab("nUMI") + ylab("nGene") + scale_x_log10() + scale_y_log10()
	print(scatter)
	dev.off()
}

if(snakemake@params[["scatter_nCount_nFeature_colHb"]]){
        jpeg(sprintf("%s/%s_ScattPlot_nCount_nFeature_colHb.jpg", outdir, label), width = 8, height = 6, units="in", res=300)
	p=ggplot(data=seurat_obj@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=hbgeneCountRatio)) + geom_point() + xlab("nUMI") + ylab("nGene") +
		scale_colour_gradient(low = "gray90", high = "black") + stat_smooth() + scale_x_log10() + scale_y_log10() + theme(legend.position="bottom")
        print(p)
        dev.off()
}

if(snakemake@params[["scatter_nCount_nFeature_perSample"]]){
	jpeg(sprintf("%s/%s_ScattPlot_nCount_nFeature_colMito_perSample.jpg", outdir, label), width = 30, height = 10, units="in", res=300) 
	p=ggplot(data=seurat_obj@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=mitoCountRatio)) + geom_point(size=0.1) + xlab("nUMI") + ylab("nGene") +
		scale_colour_gradient(low = "gray90", high = "black") + stat_smooth() + scale_x_log10() + scale_y_log10() + theme(legend.position="bottom") +  
		facet_wrap(~orig.ident,ncol=9,scales="free")
	print(p)
	dev.off()
}
