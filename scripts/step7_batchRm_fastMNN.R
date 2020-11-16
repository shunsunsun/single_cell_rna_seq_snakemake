options(stringsAsFactors = F)
suppressPackageStartupMessages({
	library(sctransform)
	library(Seurat)
	library(SeuratWrappers)
})

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
nfeatures <- as.numeric(snakemake@params[["nfeatures"]])
set.seed(1129)

se <- readRDS(file=infile)
if(!grepl("sct", snakemake@params[["norm4fastMNN"]])){
	se <- NormalizeData(se)
	se <- FindVariableFeatures(se,nfeatures=nfeatures)
}
se@meta.data$orig.ident=sapply(strsplit(se@meta.data$orig.ident,"-"),"[",1)
#se=RunFastMNN(object.list = SplitObject(se,split.by = "orig.ident"),features=nfeatures,
#	auto.merge=TRUE,BPPARAM=BiocParallel::MulticoreParam(min(as.numeric(snakemake@threads),parallel::detectCores())))
se=RunFastMNN(object.list = SplitObject(se,split.by = "orig.ident"),features=nfeatures,auto.merge=TRUE)
saveRDS(se, file=outfile)
