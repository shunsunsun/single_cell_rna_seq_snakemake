suppressPackageStartupMessages({
	library(Seurat)
	library(ggplot2)
	library(patchwork)
})

infile=snakemake@input[[1]] #e.g., data/GEJ_EP_pipeCompflted_clustered.rds
outfile=snakemake@output[[1]] #e.g., /plot/clustering/GEJ_EP_clustCompositionPlot.jpg
dimReduct=snakemake@params[[1]] #glmpca10 for GEJ, glmpca20 for ESCC
outdir=gsub("clustCompositionPlot.jpg","samples",outfile)
if(!dir.exists(outdir)){
	dir.create(outdir,recursive=T)
}
rm(list=ls())


se=readRDS(file=infile)   
Idents(object=se) <- se@meta.data$orig.ident
my_color_palette <- DiscretePalette(n=length(unique(levels(se))), palette = "alphabet")
#Use palette='polychrome' for a large number of clusters
plots <- lapply(
  X = levels(se),
  FUN = function(x) {
    return(DimPlot(
      se, group.by="seurat_clusters",
      reduction = paste0(dimReduct,"_tsne"),cols= my_color_palette,
      cells.highlight = CellsByIdentities(se, idents = x)
    ))
  }
)
p1=DimPlot(se, group.by="seurat_clusters", label=T, label.size=2, cols= my_color_palette, reduction=paste0(dimReduct,"_tsne"))
for(i in 1:length(plots)){
    p2=plots[[i]]+NoLegend()
    ggsave(file=paste0(outdir,levels(se)[i],"_tsne.jpg"), plot = p1+p2, width = 25, height = 10, device="jpeg",units="cm",dpi=300)
}

##Stack bar plot
my_color_palette <- DiscretePalette(n=length(unique(levels(se))), palette = "alphabet") 
#my_color_palette <- DiscretePalette(n=length(unique(levels(se))), palette = "polychrome") #for large number of clusters
p1=DimPlot(se, group.by="seurat_clusters", label=T, label.size=8, cols= my_color_palette, reduction=paste0(dimReduct,"_tsne"))
p2=DimPlot(se, group.by="seurat_clusters", label=T, label.size=8, cols= my_color_palette, reduction=paste0(dimReduct,"_umap"))
p = p1+p2+ plot_layout(guides ='collect') & theme(legend.text=element_text(size=18))
#basic plot of clusters by replicate
#ggplot(se@meta.data, aes(x=orig.ident, fill=seurat_clusters)) + geom_bar()
#plot as proportion or percentage of cluster
p3=ggplot(se@meta.data, aes(x=orig.ident, fill=seurat_clusters)) + geom_bar(position = "fill") + 
	scale_fill_manual(values = my_color_palette)+ theme(legend.position="none")
ggsave(file=outfile, plot = p/p3, width = 80, height = 40, device="jpeg",units="cm",dpi=300)
