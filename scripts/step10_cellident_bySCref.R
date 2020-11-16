suppressPackageStartupMessages({
	library(Seurat)
	library(SingleR)
	library(BiocParallel)
	library(pheatmap)
	library(ggplot2)
	library(patchwork)
	library(sankeywheel)
})

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
plotPrefix <- gsub("_sankeyPlot2.jpg", "", outfile)
SingleRres <- snakemake@params[["interimRes"]]
ref_file <- snakemake@params[["ref"]]
reducPrefix <- snakemake@params[["reducPrefix"]]
#nworkers=min(as.numeric(snakemake@threads),ifelse(grepl("ESCC_EP",basename(infile)), 10, parallel::detectCores()))
nworkers=4
bpParam=BiocParallel::MulticoreParam(nworkers)
#bpParam=ifelse(grepl("ESCC_EP",basename(infile)),BiocParallel::SerialParam(),BiocParallel::MulticoreParam(nworkers))

##-----Already done: Prepare locally accessible references-----
#see https://www.bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html#3_Using_single-cell_references
if(grepl("schcl", tolower(outfile))){
	if(file.exists(ref_file)){
		ref <- readRDS(file=ref_file)
	}else{
		load(file="resources/scHCL_expr.rda")
		library(SingleCellExperiment)
		library(scuttle)
		cell.labels <- colnames(ref.expr)
		gene.labels <- rownames(ref.expr)
		ref <- SingleCellExperiment(assays=list(counts=ref.expr),
    			colData=DataFrame(cell_names=colnames(ref.expr)),
			rowData=DataFrame(gene_names=rownames(ref.expr)),
			metadata=list(study="scHCL"))
		ref <- logNormCounts(ref)
		saveRDS(ref, file=ref_file)
	}
}

se <- readRDS(file=infile)
if(!file.exists(SingleRres)){
	test=GetAssayData(object = se[["SCT"]], slot = "data")
	#testidx=1:2000
	#test=test[,testidx]	
	pred <- SingleR(test = test, ref = ref, labels = ref$cell_names, de.method="wilcox", BPPARAM=bpParam)
	save(pred, file=SingleRres)
}else{
	load(file=SingleRres)
}

se[["scHCL_SingleRlabel"]]=pred$labels
se[["scHCL_SingleRscore"]]=pred$scores
saveRDS(se,file=infile)

#clust = df2%>%group_by(seurat_clusters)%>%summarise(clustSize=n())%>%as.data.frame()
#scHCLpred = df2%>%group_by(seurat_clusters,scHCL_label)%>%count(scHCL_label)%>%arrange(seurat_clusters,desc(n))%>%group_by(seurat_clusters)%>%top_n(1)%>%as.data.frame()
#singleRpred = df2%>%select(seurat_clusters,SingleR_cluster_label)%>%arrange(seurat_clusters)%>%as.data.frame()%>%unique()
#compare = merge(clust,merge(scHCLpred,singleRpred))
#colnames(compare)[4]="scHCLlabel_n"
#saveRDS(compare,file=outfile)

##per-label distribution of the deltas across cells; Labels with especially low deltas may warrant some additional caution in their interpretation
#deltaDist <- paste0(plotPrefix,"_deltaDist.jpg")
#if(!file.exists(deltaDist)){
#    jpeg(paste0(plotPrefix,"_deltaDist.jpg"), width = 50, height = 50, units="cm", res=300)
#    plotDeltaDistribution(pred,ncol=10)
#    dev.off()
#}

#combineScoreHeatmap <- paste0(plotPrefix,"_combineScoreHeatmap.jpg")
#if(!file.exists(combineScoreHeatmap)){
#	jpeg(combineScoreHeatmap, width=50, height=50, units="cm", res=300)
#	plotScoreHeatmap(pred, labels.use="labels", clusters=se$seurat_clusters, scores.use=0, calls.use=0)
#	dev.off()
#}

clustLabelHeatmap <- paste0(plotPrefix,"_clustLabelHeatmap1.jpg")
if(!file.exists(clustLabelHeatmap)){
    tab <- table(label=se$scHCL_label, cluster=se$seurat_clusters)
    # Adding a large pseudo-count of 10 to avoid strong color jumps with just 1 cell
    jpeg(clustLabelHeatmap, width = 20, height = 80, units="cm", res=300)
    #pheatmap::pheatmap(log10(tab+10), treeheight_row = 0, treeheight_col = 0, fontsize = 7)
    pheatmap::pheatmap(log10(tab+10), fontsize = 8)
    dev.off()
}
clustLabelHeatmap <- paste0(plotPrefix,"_clustLabelHeatmap2.jpg")
if(!file.exists(clustLabelHeatmap)){
    tab <- table(label=se$scHCL_SingleRlabel, cluster=se$seurat_clusters)
    # Adding a large pseudo-count of 10 to avoid strong color jumps with just 1 cell
    jpeg(clustLabelHeatmap, width = 20, height = 80, units="cm", res=300)
    #pheatmap::pheatmap(log10(tab+10), treeheight_row = 0, treeheight_col = 0, fontsize = 7)
    pheatmap::pheatmap(log10(tab+10), fontsize = 8)
    dev.off()
}

sankeyPlot <- paste0(plotPrefix,"_sankeyPlot1.jpg")
if(!file.exists(sankeyPlot)){
        jpeg(sankeyPlot, width=50, height=50, units="cm", res=300)
        sankeywheel(from=se@meta.data$seurat_clusters, to=se@meta.data$scHCL_label,
                weight=se@meta.data$scHCL_score, type="sankey", width="100%")
}
sankeyPlot <- paste0(plotPrefix,"_sankeyPlot2.jpg")
if(!file.exists(sankeyPlot)){
        jpeg(sankeyPlot, width=50, height=50, units="cm", res=300)
        sankeywheel(from=se@meta.data$seurat_clusters, to=se@meta.data$scHCL_SingleRlabel,
                weight=se@meta.data$scHCL_SingleRscore, type="sankey", width="100%")
}
