options(stringsAsFactors = F)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scry))

infile <- snakemake@input[[1]]
outPrefix <- snakemake@params[["prefix"]]
nHGV <- as.numeric(snakemake@params[["nhgv"]])
plotDir <- snakemake@params[["plot"]]

se <- readRDS(file=infile)
for(assay in c("SCT","RNA")){
	sce=as.SingleCellExperiment(se, assay=assay)
	sce<-devianceFeatureSelection(sce, assay="counts", sorted=TRUE)
	jpeg(sprintf("%s/%s_%scount_genedevRank.jpg",plotDir,gsub("_sctNorm.rds","",basename(infile)),assay),width=10,height=5,unit="in",res=300)		
	p=plot(rowData(sce)$binomial_deviance, type="l", xlab="ranked genes",
     		ylab="binomial deviance", main="Feature Selection with Deviance")
	print(p)
	dev.off()
	#tmp=rowData(sce)[["binomial_deviance"]]
	#hgv=names(tmp[order(tmp, decreasing=T)[1:min(nHGV,length(tmp))]])
	sce2<-GLMPCA(sce[1:nHGV,], assay="counts", 2)
	fit<-metadata(sce2)$glmpca
	save(sce2,fit,file=paste0(outPrefix,"_glmpca_",assay,"count.rda")
	#pd<-cbind(as.data.frame(colData(sce2)), fit$factors)
	#ggplot(pd, aes(x=dim1, y=dim2, colour=phenoid)) + geom_point(size=.8) +
	#	ggtitle("GLM-PCA applied to high deviance genes")
}
