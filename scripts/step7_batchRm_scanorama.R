options(stringsAsFactors = F)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(reticulate))
scanorama <- import('scanorama')

infile <- snakemake@input[[1]]
if(grepl("sct", snakemake@params[["norm4harmony"]])){
	targetAssay <- "SCT"
}else{
	targetAssay <- "RNA"
}
outfile <- snakemake@output[["rds"]]
theta <- as.numeric(snakemake@params[["theta"]])
nclust = as.numeric(snakemake@params[["nclust"]])
max_it_clust = as.numeric(snakemake@params[["max_it_clust"]])


# https://github.com/brianhie/scanorama/issues/38#issuecomment-551446738

se <- readRDS(file=infile)
se <- RunPCA(object=se, verbose=T)  ##https://satijalab.org/seurat/v3.2/sctransform_vignette.html
se <- RunHarmony(object=se, assay.use = targetAssay, reduction = "pca", dims.use = 1:50, group.by.vars = "orig.ident",
	plot_convergence = TRUE, theta=theta, nclust=nclust, max.iter.cluster=max_it_clust)

saveRDS(se, file=outfile)
