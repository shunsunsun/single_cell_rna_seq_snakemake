options(stringsAsFactors = F)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))

indir <- snakemake@params[["indir"]]
label <- snakemake@params[["infile"]]
outfile <- snakemake@output[["rds"]]

nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
#cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 25*1024^3)
set.seed(1129)

EP=scan(file=snakemake@input[["EP_list"]],what="character")
IM=scan(file=snakemake@input[["IM_list"]],what="character")
samples = union(EP,IM)
normed_seurat <- future_lapply(X = samples, FUN = function(s){
        readRDS(file=paste0(indir,"/",s,".",label,".rds"))
},=randomseed)
merged=merge(normed_seurat[[1]],y=normed_seurat[-1],add.cell.ids=samples)
saveRDS(normed_seurat, file=outfile)
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
