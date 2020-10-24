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

## Read input
se <- readRDS(infile)

## Load cell cycle markers
load(snakemake@params[["ccgenes"]])

plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 25*1024^3) #25G
set.seed(1129)

se <- NormalizeData(se)
se <- FindVariableFeatures(se, selection.method = "vst") 

tryCatch({
	tmp <- CellCycleScoring(se, g2m.features=g2m_genes, s.features=s_genes)
	tmp$cc_difference <- tmp$S.Score - tmp$G2M.Score
	se <- ScaleData(tmp, features=rownames(tmp), vars.to.regress = c(regVars, 'cc_difference'))
}, error=function(cond){  ##in case cell cycle scoring has any problem
	cat(sprintf("No cell cycle scoring due to %s",cond))
        se <- ScaleData(se, features=rownames(se), vars.to.regress = regVars)
})
saveRDS(se, file=outfile)
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
