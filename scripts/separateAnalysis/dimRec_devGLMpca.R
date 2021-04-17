suppressPackageStartupMessages({
	library(sctransform)
	library(Seurat)
	library(SeuratWrappers)
	library(glmpca)
})

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
outfile <- args[2]
nHVG <- 3000
ndims <- c(10,20) 
set.seed(1129)

se <- readRDS(file=infile)
DefaultAssay(se) <- "RNA"
se <- DietSeurat(se, assays="RNA")

m <- GetAssayData(se, slot = "counts", assay = "RNA")
devs <- scry::devianceFeatureSelection(m)
dev_ranked_genes <- rownames(se)[order(devs, decreasing = TRUE)]
for(n in ndims){
	#Sparse matrices are coerced to dense matrice for  minibatch='none'; If this exhausts memory, consider setting minibatch to 'stochastic' or 'memoized'
	se <- RunGLMPCA(se, features = head(dev_ranked_genes,n=nHVG), L = n, minibatch='stochastic',
		reduction.name=paste0('glmpca',n), reduction.key=paste0("GLMPC",n,"_"))
}
saveRDS(se,file=outfile)

#plotfile=sprintf("%s/%s_devrank.jpg",plotDir,gsub(".rds","",basename(outfile)))
#devs = devs[order(devs,decreasing=TRUE)]
#jpeg(plotfile, width = 8, height = 6, units="cm", res=300)
#plot(unname(devs), type="l",xlab="ranked genes",ylab="binomial deviance")
#abline(v=nHVG, lty=2, col="red")
#dev.off()
