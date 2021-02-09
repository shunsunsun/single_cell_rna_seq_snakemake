#source /appsnew/source/R-4.0.2share.sh
#R.Version()
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(ggplot2))

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
outPrefix <- args[2]
reduc <- args[3] #integrated (CCA based clustering) or SCT (harmony based clustering)
randSeed <- 1129L
set.seed(randSeed)

se <- readRDS(file=infile)

suppressPackageStartupMessages(library(clustree))
p <- clustree(df,prefix=paste0(reduc,"_snn_res."),edge_width=0.8,node_alpha = 0.8) + 
	scale_color_brewer(palette = "Set1") + 
	theme(legend.position = "bottom") + 
	scale_edge_color_continuous(low = "grey80", high = "red")
ggsave(filename=paste0(outPrefix,"_clustree.pdf"), plot = p, width=15, height=10, units="in", dpi=300)

# Feature plot with QC metrics
metrics=c("log10_total_counts","log10_total_features","mitoCountRatio","cc_difference","featcount_ratio")
p1 <- FeaturePlot(se, reduction = "umap", features = metrics, ncol=2, pt.size = .4, order = TRUE, min.cutoff = 'q10')
ggsave(filename=paste0(outPrefix,"_qc.pdf"), plot = p1, width = 12, height = 15, units="in",dpi=300)

## UMAP grouped and colored by cell cycle phases
p2 <- DimPlot(se, reduction = "umap", group.by='Phase', split.by="Phase")
ggsave(filename=paste0(outPrefix,"_cc.pdf"), plot = p2, width = 12, height = 4, units="in",dpi=300)
