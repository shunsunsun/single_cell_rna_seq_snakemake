suppressPackageStartupMessages({
	library(Seurat)
	library(ggplot2)
	library(patchwork)
})

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
plotPrefix <- paste0(snakemake@params[["outPlot"]], '/', gsub("_MarkerPred.txt", "", basename(outfile)))
reducPrefix <- snakemake@params[["reducPrefix"]]

marker_df=read.table(file="resources/markers.tsv",header=T)
marker_list <- list()
for(r in 1:nrow(marker_df)){
	marker_list[[marker_df$celltype[r]]]=unname(unlist(strsplit(marker_df$genes[r],',')))
}

se <- readRDS(file=infile)
allgenes=rownames(se)
for(i in 1:length(marker_list)){
	celltype=names(marker_list)[i]
	markergenes=marker_list[[celltype]]
	print(paste(intersect(markergenes,allgenes),collapse=","))
	dot_plot <- paste0(plotPrefix,"_",celltype,"_dotPlot.pdf")
	#if(!file.exists(dot_plot)){
		p <- DotPlot(se, features=markergenes,group.by='seurat_clusters')+coord_flip()
		ggsave(filename=dot_plot, plot=p, width=15, height=3*length(markergenes),unit="in",dpi=300)
	#}
	feat_plot <- paste0(plotPrefix,"_",celltype,"_featurePlot.pdf")
	#if(!file.exists(feat_plot)){		
		p <- FeaturePlot(se,reduction=paste0(reducPrefix,"_umap"),features=markergenes,order=T,min.cutoff='q10',repel=T,ncol=min(length(markergenes),2))
		ggsave(filename=feat_plot, plot=p, width=10, height=ceiling(length(markergenes)/2)*4,unit="in",dpi=300)
	#}
	vln_plot <- paste0(plotPrefix,"_",celltype,"_vlnPlot.pdf")
	#if(!file.exists(vln_plot)){
		p <- VlnPlot(se,features=markergenes,ncol=min(length(markergenes),2))
		ggsave(filename=vln_plot, plot=p, width=10, height=ceiling(length(markergenes)/2)*4,unit="in",dpi=300)
	#}
}
write("Done",file=outfile)
