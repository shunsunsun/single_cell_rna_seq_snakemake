options(stringsAsFactors = F)
library(Seurat)

samples=scan(file=snakemake@input[["list"]],what="character")
sceList=lapply(samples, function(s){
	readRDS(file=paste0(snakemake@params[[1]],"/",s,".seuratobj.rds"))
})

#Because the same cell IDs can be used for different samples, we add a sample-specific prefix to each of our cell IDs using the add.cell.id argument
merged=merge(sceList[[1]],y=sceList[-1],add.cell.ids=samples,project="gej_escc_sc")
saveRDS(merged,file=snakemake@output[[1]])
