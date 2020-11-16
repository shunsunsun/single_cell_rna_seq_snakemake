suppressPackageStartupMessages({
	library(Seurat)
	library(SingleR)
	library(celldex)
	library(BiocParallel)
	library(pheatmap)
	library(ggplot2)
	library(patchwork)
})

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
plotPrefix <- gsub("_SingleRlabelTSNEUMAP.jpg", "", outfile)
singleRres <- snakemake@params[["interimRes"]]
ref <- snakemake@params[["ref"]]
label <- snakemake@params[["label"]]  #label.main or label.fine; the latter represents the highest resolution of annotation available for a reference set
refdir <- snakemake@params[["refdir"]]
reducPrefix <- snakemake@params[["reducPrefix"]]
nworkers=min(as.numeric(snakemake@threads),ifelse(grepl("ESCC_EP",basename(infile)), 10, parallel::detectCores()))
bpParam=BiocParallel::MulticoreParam(nworkers)
#bpParam=ifelse(grepl("ESCC_EP",basename(infile)),BiocParallel::SerialParam(),BiocParallel::MulticoreParam(nworkers))


##-----Already done: Prepare locally accessible references-----
#Run the following code at the login node to download and process reference data as our computing nodes have no internet access
#The following references in celldex are all for humans
#See https://bioconductor.org/packages/devel/data/experiment/vignettes/celldex/inst/doc/userguide.html
if(!dir.exists(refdir)){
	dir.create(refdir,recursive=T)
	hpca <- HumanPrimaryCellAtlasData() #include cells/cell lines treated or collected from pathogenic conditions
	saveRDS(hpca, file=paste0(refdir,'/', 'HumanPrimaryCellAtlasData.rds'))
	be <- BlueprintEncodeData() #suitable for quick and interpretable annotation of stroma and immune cells
	saveRDS(be, file=paste0(refdir,'/', 'BlueprintEncodeData.rds'))
	dice <- DatabaseImmuneCellExpressionData() #rich CD4+ T cell subpopulations but lacking dendritic cells and only having one B cell subset (naive)
	saveRDS(dice, file=paste0(refdir,'/', 'DatabaseImmuneCellExpressionData.rds'))
	nh <- NovershternHematopoieticData() #greatest resolution for myeloid and progenitor cells; fewer T cell subsets than other immune references but many more NK, erythroid, and granulocytic subsets
	saveRDS(nh, file=paste0(refdir,'/', 'NovershternHematopoieticData.rds'))
	mi <- MonacoImmuneData() #expansive B and T cell subsets
	saveRDS(mi, file=paste0(refdir,'/', 'MonacoImmuneData.rds'))

	##Utilize cell ontology
	library(ontoProc)
	bfc <- BiocFileCache::BiocFileCache(ask=FALSE)
	path <- BiocFileCache::bfcrpath(bfc, "http://purl.obolibrary.org/obo/cl.obo")
	cl <- get_ontology(path, extract_tags="everything")
	parents <- cl$parents
	self <- rep(names(parents), lengths(parents))
	library(igraph)
	g <- make_graph(rbind(unlist(parents), self))
	save(cl, g, file=paste0(refdir,'/', 'cell_ontology_index.rda'))
}

##--------Get reference---------

ref_list <- list()
label_list <- list()
#co_label_list <- list()
if(grepl("hpca",ref)){
	ref_list[["hpca"]] <- readRDS(file=paste0(refdir,'/', 'HumanPrimaryCellAtlasData.rds'))
	label_list[["hpca"]] <- ref_list[["hpca"]][[label]]
#	co_label_list[["hpca"]] <- ref_list[["hpca"]][["label.ont"]]
}
if(grepl("be",ref)){
	ref_list[["be"]] <- readRDS(file=paste0(refdir,'/', 'BlueprintEncodeData.rds'))
	label_list[["be"]] <- ref_list[["be"]][[label]]
#	co_label_list[["be"]] <- ref_list[["be"]][["label.ont"]]
}
if(grepl("dice",ref)){
	ref_list[["dice"]] <- readRDS(file=paste0(refdir,'/', 'DatabaseImmuneCellExpressionData.rds'))
	label_list[["dice"]] <- ref_list[["dice"]][[label]]
#	co_label_list[["dice"]] <- ref_list[["dice"]][["label.ont"]]
}
if(grepl("nh",ref)){
	ref_list[["nh"]] <- readRDS(file=paste0(refdir,'/', 'NovershternHematopoieticData.rds'))
	label_list[["nh"]] <- ref_list[["nh"]][[label]]
#	co_label_list[["nh"]] <- ref_list[["nh"]][["label.ont"]]
}
if(grepl("mi",ref)){
	ref_list[["mi"]] <- readRDS(file=paste0(refdir,'/', 'MonacoImmuneData.rds'))
	label_list[["mi"]] <- ref_list[["mi"]][[label]]
#	co_label_list[["mi"]] <- ref_list[["mi"]][["label.ont"]]
}
##!Alread done: Consistency in the main labels of these references
# for(k in 1:(length(ref_list)-1)){
    # for(j in (k+1):length(ref_list)){
        # name1 <- names(ref_list)[k]
        # fontsize = ifelse(name1=="hpca", 7, 12)
        # name2 <- names(ref_list)[j]
        # filename <- paste0(refdir,'/RefLabelCompare_',name1,'_',name2,'.jpg')
        # if(!file.exists(filename)){
            # matched <- matchReferences(ref_list[[name1]], ref_list[[name2]], ref_list[[name1]]$label.main, ref_list[[name2]]$label.main)
            # jpeg(filename, width = 10, height = 10, units="cm", res=300)
            # pheatmap::pheatmap(matched, col=viridis::plasma(100), treeheight_row = 0, treeheight_col = 0, fontsize = fontsize)
            # dev.off()
        # }
    # }
# }


if(!file.exists(singleRres)){
	##-------Get test data------------
	se <- readRDS(file=infile)
	test=GetAssayData(object = se[["SCT"]], slot = "data")
	#testidx=1:2000
	#test=test[,testidx]
	#finalClust=finalClust[testidx]


	##------SingleR identification--------

	pred.classic <- SingleR(test = test, ref = ref_list, labels = label_list, recompute=TRUE, BPPARAM=bpParam)
	print("finish")
	#pred.colabel <- SingleR(test = test, ref = ref_list, labels = co_label_list, recompute=TRUE, BPPARAM=bpParam)

	##Alternative statistical testing for marker gene detection: Wilcoxon ranked sum test and Weltch t test, both considering the variance of gene expression across cells
	##These tests are slower but more appropriate for single-cell reference data compared to the default marker detection algorithm, which may fail for low-coverage data where the median for each label is often zero.
	#pred.wilcox <- SingleR(test = test, ref = ref_list, labels = label_list, recompute=TRUE, de.method="wilcox", BPPARAM=bpParam)

	#Annotation can be performed on the cluster-level profiles rather than on the single-cell level. it directly returns the likely cell type identity of each cluster
	#This approach assumes that each cluster in the test dataset corresponds to exactly one reference label, which is likely not true
	clust.pred <- SingleR(test = test, ref = ref_list, labels = label_list, recompute=TRUE, clusters=finalClust, BPPARAM=bpParam)
	print("finish")
	#clust.pred.colabel <- SingleR(test = test, ref = ref_list, labels = co_label_list, recompute=TRUE, clusters=finalClust, BPPARAM=bpParam)
	save(pred.classic, clust.pred, file=singleRres)
}else{
	load(file=singleRres)
	if(file.exists(gsub(".rds","_SingleRann.rds",infile))){
		se=readRDS(file=gsub(".rds","_SingleRann.rds",infile))
	}else{
		se=readRDS(file=infile)
	}
}

##per-label distribution of the deltas across cells; Labels with especially low deltas may warrant some additional caution in their interpretation
if(!file.exists(paste0(plotPrefix,"_deltaDist1.jpg"))){
    plist <- plotDeltaDistribution(pred.classic,grid.vars=NULL,ncol=10)
    for(i in 1:length(plist)){
        jpeg(paste0(plotPrefix,"_deltaDist",i,".jpg"), width = 20, height = 80, units="cm", res=300)
        print(plist[[i]])
        dev.off()
    }
}

combineScoreHeatmap <- paste0(plotPrefix,"_SingleRcombineScoreHeatmap.jpg")
if(!file.exists(combineScoreHeatmap)){
	jpeg(combineScoreHeatmap, width=50, height=50, units="cm", res=300)
	plotScoreHeatmap(pred.classic, labels.use="pruned.labels", clusters=se$seurat_clusters, scores.use=0, calls.use=0)
	dev.off()
}

clustLabelHeatmap <- paste0(plotPrefix,"_SingleRclustLabelHeatmap.jpg")
if(!file.exists(clustLabelHeatmap)){
    tab <- table(label=pred.classic$labels, cluster=se$seurat_clusters)
    # Adding a large pseudo-count of 10 to avoid strong color jumps with just 1 cell
    jpeg(clustLabelHeatmap, width = 20, height = 80, units="cm", res=300)
    #pheatmap::pheatmap(log10(tab+10), treeheight_row = 0, treeheight_col = 0, fontsize = 7)
    pheatmap::pheatmap(log10(tab+10), fontsize = 8)
    dev.off()
}

if(!file.exists(gsub(".rds","_SingleRann.rds",infile))){
	se[["SingleR_label"]]=pred.classic$labels
	scores = as.data.frame(pred.classic$scores)
	se[["SingleR_score"]]=apply(scores, 1, function(x){max(x,na.rm=T)})
	se[["SingleR_cluster_label"]] <- clust.pred$labels[match(se@meta.data$seurat_clusters, rownames(clust.pred))]
	saveRDS(se,file=gsub(".rds","_ann.rds",infile))
}

library(sankeywheel)
sankeyPlot <- paste0(plotPrefix,"_SingleRsankeyPlot.jpg")
if(!file.exists(sankeyPlot)){
        jpeg(sankeyPlot, width=50, height=50, units="cm", res=300)
        sankeywheel(from=se@meta.data$seurat_clusters, to=se@meta.data$SingleR_label, 
		weight=se@meta.data$SingleR_score, type="sankey", width="100%")
}

tsne_umap_label <- outfile
if(!file.exists(tsne_umap_label)){
    p1 = DimPlot(se, group.by="SingleR_cluster_label", label=T, label.size=2, reduction=paste0(reducPrefix, "_tsne"))
    p2 = DimPlot(se, group.by="SingleR_cluster_label", label=T, label.size=2, reduction=paste0(reducPrefix, "_umap"))
    p = p1+p2+ plot_layout(guides ='collect') & theme(legend.text=element_text(size=8))
    ggsave(file=tsne_umap_label, plot = p, width = 25, height = 10, device="jpeg",units="cm",dpi=300)
}
