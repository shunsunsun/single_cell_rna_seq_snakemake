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
	dot_plot <- paste0(plotPrefix,"_",celltype,"_dotPlot.jpg")
	#if(!file.exists(dot_plot)){
		jpeg(file=dot_plot,width=30,height=4*length(markergenes),unit="cm",res=300)
		DotPlot(se, features=markergenes,group.by='seurat_clusters')+coord_flip()
		dev.off()
	#}
	feat_plot <- paste0(plotPrefix,"_",celltype,"_featurePlot.jpg")
	#if(!file.exists(feat_plot)){
		jpeg(file=feat_plot,width=50,height=ceiling(length(markergenes)/2)*20,unit="cm",res=300)
		FeaturePlot(se,reduction=paste0(reducPrefix,"_umap"),features=markergenes,order=T,min.cutoff='q10',repel=T,ncol=min(length(markergenes),2))
		dev.off()	
	#}
	vln_plot <- paste0(plotPrefix,"_",celltype,"_vlnPlot.jpg")
	#if(!file.exists(vln_plot)){
		jpeg(file=vln_plot,width=50,height=ceiling(length(markergenes)/2)*20,unit="cm",res=300)
		VlnPlot(se,features=markergenes,ncol=min(length(markergenes),2))
		dev.off()
	#}
}
write("Done",file=outfile)
