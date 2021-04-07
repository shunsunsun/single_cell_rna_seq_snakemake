suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(future))

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
outfile <- args[2]
outdir <- dirname(outfile)
if(!dir.exists(outdir)){
        dir.create(outdir,recursive=T)
}

regVars=c('mitoCountRatio', 'nFeature_RNA', 'nCount_RNA', 'percent_mito', 'percent.mt')

nworker <- min(as.numeric(args[3]),length(availableWorkers()))
#cat(sprintf("Use %d workders\n",nworker))

## Read input
flted_se <- readRDS(infile)
regVars=intersect(regVars,colnames(flted_se[[]]))

## Load cell cycle markers
load(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/resources/cycle.rda")

plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 20*1024^3) #20G
set.seed(1129L)

# https://github.com/immunogenomics/harmony/issues/41#issuecomment-633885490
# if samples from different technical platforms, try https://github.com/immunogenomics/harmony/issues/41#issuecomment-642862186
tryCatch({
	tmp <- SCTransform(flted_se, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = regVars)
        tmp <- CellCycleScoring(tmp, s.features = s_genes, g2m.features = g2m_genes, assay = 'SCT', set.ident = TRUE)
	tmp$cc_difference <- tmp$S.Score - tmp$G2M.Score
        normed_se <- SCTransform(tmp, assay = 'RNA', new.assay.name = 'SCT', variable.features.n=5000, vars.to.regress = c(regVars, 'cc_difference'))
}, error=function(cond){
	cat(sprintf("No cell cycle scoring due to %s",cond))
        normed_se <- SCTransform(flted_se, variable.features.n=5000, vars.to.regress = regVars)
})
saveRDS(normed_se, file=outfile)
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
