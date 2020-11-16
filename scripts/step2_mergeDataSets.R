options(stringsAsFactors = F)
suppressPackageStartupMessages(library(Seurat))

projectName=snakemake@params[[1]]
filePaths=snakemake@input[[1]]
outfile=snakemake@output[[1]]

datasets=scan(file=filePaths,what="character")
sceList=lapply(datasets, function(d){
	readRDS(file=d)
})

merged=merge(sceList[[1]],y=sceList[-1],project=projectName)
saveRDS(merged,file=outfile)
