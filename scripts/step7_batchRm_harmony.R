options(stringsAsFactors = F)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(harmony))

infile <- snakemake@input[[1]]
#if(grepl("sct", snakemake@params[["norm4harmony"]])){
	targetAssay <- "SCT"
#}else{
#	targetAssay <- "RNA"
#}
outfile <- snakemake@output[[1]]
#theta <- as.numeric(snakemake@params[["theta"]])
#nclust = as.numeric(snakemake@params[["nclust"]])
#max_it_clust = as.numeric(snakemake@params[["max_it_clust"]])


# https://github.com/immunogenomics/harmony/issues/41#issuecomment-633885490
# if samples from different techical platforms, try https://github.com/immunogenomics/harmony/issues/41#issuecomment-642862186

se <- readRDS(file=infile)
se@meta.data$orig.ident=sapply(strsplit(se@meta.data$orig.ident,"-"),"[",1)
se <- RunPCA(object=se, verbose=T)  ##https://satijalab.org/seurat/v3.2/sctransform_vignette.html
#se <- RunHarmony(object=se, assay.use = targetAssay, reduction = "pca", dims.use = 1:50, group.by.vars = "orig.ident",
#	plot_convergence = FALSE, theta=theta, nclust=nclust, max.iter.cluster=max_it_clust)
se <- RunHarmony(object=se, assay.use = targetAssay, reduction = "pca", dims.use = 1:50, group.by.vars = "orig.ident", plot_convergence = FALSE)
saveRDS(se, file=outfile)
