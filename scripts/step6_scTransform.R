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
flted_se <- readRDS(infile)

## Load cell cycle markers
load(snakemake@params[["ccgenes"]])

plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 20*1024^3) #20G
set.seed(1129)

# https://github.com/immunogenomics/harmony/issues/41#issuecomment-633885490
# if samples from different techical platforms, try https://github.com/immunogenomics/harmony/issues/41#issuecomment-642862186
if(!snakemake@params[["sctPreNorm"]]){
	tryCatch({
               	tmp <- NormalizeData(flted_se)
		tmp <- CellCycleScoring(tmp, g2m.features=g2m_genes, s.features=s_genes)
		tmp$cc_difference <- tmp$S.Score - tmp$G2M.Score
		normed_se <- SCTransform(tmp, variable.features.n=5000, vars.to.regress = c(regVars, 'cc_difference'))
	}, error=function(cond){  ##in case cell cycle scoring has any problem
               	#if(grepl("Insufficient data values to produce", cond, fixed=TRUE)){
		cat(sprintf("No cell cycle scoring due to %s",cond))
                normed_se <- SCTransform(flted_se, variable.features.n=5000, vars.to.regress = regVars)
                #}
        })
}else{
        tryCatch({
               	tmp <- SCTransform(flted_se, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = regVars)
               	tmp <- CellCycleScoring(tmp, s.features = s_genes, g2m.features = g2m_genes, assay = 'SCT', set.ident = TRUE)
		tmp$cc_difference <- tmp$S.Score - tmp$G2M.Score
                normed_se <- SCTransform(tmp, assay = 'RNA', new.assay.name = 'SCT', variable.features.n=5000, vars.to.regress = c(regVars, 'cc_difference'))
       	}, error=function(cond){
		cat(sprintf("No cell cycle scoring due to %s",cond))
                normed_se <- SCTransform(flted_se, variable.features.n=5000, vars.to.regress = regVars)
        })
}
saveRDS(normed_se, file=outfile)
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
