suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratWrappers))
suppressPackageStartupMessages(library(glmpca))
suppressPackageStartupMessages(library(scry))

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
nHVG <- as.numeric(snakemake@params[["nhvg"]])
plotDir <- snakemake@params[["plot"]]
ndim <- as.numeric(snakemake@params[["ndim"]])
set.seed(1129)

se <- readRDS(file=infile)
m <- GetAssayData(se, slot = "counts", assay = "RNA")
devs <- scry::devianceFeatureSelection(m)
dev_ranked_genes <- rownames(se)[order(devs, decreasing = TRUE)]
for(n in c(2000,nHVG)){
	topdev <- head(dev_ranked_genes, n)
	#warning:Sparse matrices are not supported for minibatch='none'. Coercing to dense matrix. If this exhausts memory, consider setting minibatch to 'stochastic' or 'memoized'.
	se <- RunGLMPCA(se, features = topdev, L = ndim, minibatch='stochastic', reduction.name=paste0('glmpca_dev',n), reduction.key=paste0("dev",n,"GLMPC_"))
}
saveRDS(se,file=outfile)

plotfile=sprintf("%s/%s_devrank.jpg",plotDir,gsub(".rds","",basename(outfile)))
devs = devs[order(devs,decreasing=TRUE)]
jpeg(plotfile, width = 10, height = 8, units="in", res=300)
plot(unname(devs), type="l",xlab="ranked genes",ylab="binomial deviance",
	main=paste0("Feature selection with deviance (latent space ndim = ", ndim, ")"))
abline(v=2000, lty=2, col="red")
abline(v=nHVG, lty=2, col="red")
dev.off()
