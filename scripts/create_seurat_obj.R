options(stringsAsFactors = F)
library(Seurat)

data=Read10X(data.dir=snakemake@params[["dir"]])
sample=snakemake@params[["id"]]
obj=CreateSeuratObject(counts = data, project = sample, min.cells=3, min.features=200)
saveRDS(obj,file=snakemake@output[[1]])
