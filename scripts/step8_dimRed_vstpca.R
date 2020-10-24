options(stringsAsFactors = F)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scater))

infile <- snakemake@input[[1]]
outfile <- snakemake@output[["rds"]]
report <- snakemake@output[["optimPC"]]
plotDir <- snakemake@params[["plot"]]
nHVG <- as.numeric(snakemake@params[["nhvg"]])
assay <- ifelse(grepl("sct",infile),"SCT", "RNA")

se <- readRDS(file=infile)
hgv=VariableFeatures(se, assay=assay)
hgv=hgv[1:min(nHVG,length(hgv))]

se <- RunPCA(object=se, features=hgv, reduction.name="pca_4000hvg", reduction.key = "4000hvgPC_")
se <- RunPCA(object=se, features=hgv[1:min(2000,length(hgv))], reduction.name="pca_2000hvg", reduction.key = "2000hvgPC_")
saveRDS(se, file=outfile)

#########Identify significant PCs#############

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
pca_plot=sprintf("%s/%s_elbow.jpg",plotDir,gsub(".rds","",basename(outfile)))
p1=ElbowPlot(se,reduction="pca_4000hvg",ndim=50)
p2=ElbowPlot(se,reduction="pca_2000hvg",ndim=50)
ggsave(file=pca_plot, plot = p1+p2, width = 12, height = 5, device="jpeg",units="in",dpi=300)

##Elbow quantitative metric: from pipeComp/R/misc.R
#' Identifies the point farthest from a line passing through by the first and last points. Used for automatization of the elbow method.
#' @param y Monotonically inscreasing or decreasing values
#' @param x Optional x coordinates corresponding to `y` (defaults to seq)
#' @return The value of `x` farthest from the diagonal.
farthestPoint <- function(y, x=NULL){
  if(is.null(x)) x <- seq_len(length(y))
  d <- apply( cbind(x,y), 1, 
              a=c(1,y[1]), b=c(length(y),rev(y)[1]), 
              FUN=function(y, a, b){
                v1 <- a-b
                v2 <- y-a
                abs(det(cbind(v1,v2)))/sqrt(sum(v1*v1))
              })
  order(d,decreasing=TRUE)[1]
}

##Elbow quantiative metric: from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
#We can calculate where the principal components start to elbow by taking the smaller value of:
#The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
#The point where the percent change in variation between the consecutive PCs is less than 0.1%.

hbc <- function(sdv){
	# Determine percent of variation associated with each PC
	pct <- sdv / sum(sdv) * 100
	# Calculate cumulative percents for each PC
	cumu <- cumsum(pct)
	# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
	co1 <- which(cumu > 90 & pct < 5)[1]
	# Determine the difference between variation of PC and subsequent PC. The last point where change of % of variation is more than 0.1%.
	co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
	# Minimum of the two calculation
	min(co1, co2)
}

suppressPackageStartupMessages(library(intrinsicDimension))
res=NULL
for(n in c(2000,nHVG)){
	label=paste0("pca_",n,"hvg")
	x <- Embeddings(se[[label]])
	sdv <- Stdev(se, label)
	maxLikGlobal20=maxLikGlobalDimEst(x, k=20, unbiased=TRUE)
	maxLikGlobal10=maxLikGlobalDimEst(x, k=10, unbiased=TRUE)
	elbow=farthestPoint(sdv)-1
	elbow2=hbc(sdv)
	npc=as.integer(round(max(maxLikGlobal20$dim.est,maxLikGlobal10$dim.est,elbow,elbow2)))
	res=rbind(res,data.frame(nhvg=n,maxGlob20=maxLikGlobal20$dim.est,maxGlob10=maxLikGlobal10$dim.est,elbow=elbow,elbow2=elbow2,optimPC=npc))
	rm(x,sdv,maxLikGlobal20,maxLikGlobal10,elbow,elbow2,npc)
}
write.table(res,file=report,row.names=F,col.names=T,sep="\t",quote=F)
