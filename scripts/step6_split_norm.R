suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))

input <- snakemake@input[[1]] ## filtered merged data as a seurat obj
outfile <- snakemake@output[[1]]
outdir <- snakemake@params[["dir"]]

if(!dir.exists(outdir)){
        dir.create(outdir,recursive=T)
}
outdir1 <- dirname(outfile)
if(!dir.exists(outdir1)){
        dir.create(outdir1,recursive=T)
}
suffix=gsub("Obj","",basename(outdir))

#It is technically unsound to regress out nUMI and nGene, which scTransform automatically regresses out
if(snakemake@params[["regressNum"]]){
    regVars=c('mitoCountRatio', 'nFeature_RNA', 'nCount_RNA')
}else{
    regVars=c('mitoCountRatio')
}

nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))


## Read input
flted_seurat <- readRDS(input)

## Load cell cycle markers
load(snakemake@params[["ccgenes"]])

plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 20*1024^3)

# Split seurat object to perform cell cycle scoring and SCTransform on all samples
split_seurat <- SplitObject(flted_seurat, split.by = "orig.ident")

# Normalize each sample
future_lapply(X = split_seurat, future.seed=1129, FUN = function(s) {
	sample=unique(s@meta.data$orig.ident)
	##After filtering, sample-level cell/feature number criterion may no long be met
	x=CreateSeuratObject(counts=GetAssayData(object=s, slot='counts'), meta.data=s@meta.data, min.cells=1, min.features=1)
	cat(sprintf("Processing %s with %d features and %d cells", sample, dim(x)[1], dim(x)[2]),file=outfile,append=T,sep="\n")
	if(!snakemake@params[["sctPreNorm"]]){
		tryCatch({
                	tmp <- NormalizeData(x)
                        tmp <- CellCycleScoring(tmp, g2m.features=g2m_genes, s.features=s_genes)
			tmp$cc_difference <- tmp$S.Score - tmp$G2M.Score
                        x <- SCTransform(tmp, variable.features.n=5000, vars.to.regress = c(regVars, 'cc_difference'))
		}, error=function(cond){  ##in case cell cycle scoring has any problem
                	#if(grepl("Insufficient data values to produce", cond, fixed=TRUE)){
			cat(sprintf("No cell cycle scoring due to %s",cond),file=outfile,append=T,sep="\n")
                        x <- SCTransform(x, variable.features.n=5000, vars.to.regress = regVars)
                        #}
                })
        }else{
                tryCatch({
                        tmp <- SCTransform(x, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = regVars)
                        tmp <- CellCycleScoring(tmp, s.features = s_genes, g2m.features = g2m_genes, assay = 'SCT', set.ident = TRUE)
			tmp$cc_difference <- tmp$S.Score - tmp$G2M.Score
                        x <- SCTransform(tmp, assay = 'RNA', new.assay.name = 'SCT', variable.features.n=5000, vars.to.regress = c(regVars, 'cc_difference'))
                }, error=function(cond){
			cat(sprintf("No cell cycle scoring due to %s",cond),file=outfile,append=T,sep="\n")
                        x <- SCTransform(x, variable.features.n=5000, vars.to.regress = regVars)
                })
        }
	saveRDS(x, file=paste0(outdir,"/",sample,".",suffix,".rds"))
})

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
