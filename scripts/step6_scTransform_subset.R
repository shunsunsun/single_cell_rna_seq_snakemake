#https://github.com/satijalab/seurat/issues/2014#issuecomment-629358390

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(future))

infile <- snakemake@input[[1]] ## filtered merged data as a seurat obj
outfile <- snakemake@output[[1]]
outdir <- dirname(outfile)
if(!dir.exists(outdir)){
        dir.create(outdir,recursive=T)
}

#Although scTransform automatically regresses out nUMI, at this point (after filteration), the distributions of nFeature_RNA and nCount_RNA may have been changed
#So I'm not sure whether nCount_RNA and nFeature_RNA should also be included in the variables to be regressed out. The official vignete does not but I see several people do so...
if(snakemake@params[["regressNum"]]){
    regVars=c('mitoCountRatio', 'nFeature_RNA', 'nCount_RNA')
}else{
    regVars=c('mitoCountRatio')
}

nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))

se <- readRDS(infile)

plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 20*1024^3) #20G
set.seed(1129)

se <- SCTransform(se, variable.features.n=3000, vars.to.regress = c(regVars, 'cc_difference'))
saveRDS(se, file=outfile)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
