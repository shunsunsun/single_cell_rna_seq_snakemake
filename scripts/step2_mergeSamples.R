options(stringsAsFactors = F)
suppressPackageStartupMessages(library(Seurat))

projectName=snakemake@params[["project"]]
suffix=snakemake@params[["suffix"]] #e.g., ".seuratobj.rds"
indir=snakemake@params[["indir"]]
outfile=snakemake@output[[1]]

samples=scan(file=snakemake@input[["list"]],what="character")
sceList=lapply(samples, function(s){
	readRDS(file=paste0(indir,"/",s,suffix))
})

#Because the same cell IDs can be used for different samples, we add a sample-specific prefix to each of our cell IDs using the add.cell.id argument
merged=merge(sceList[[1]],y=sceList[-1],add.cell.ids=samples,project=projectName)
saveRDS(merged,file=outfile)
