##Ref: https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
##Ref: https://satijalab.org/seurat/articles/sctransform_vignette.html

options(stringsAsFactors = F)
suppressPackageStartupMessages({
        library(sctransform)
        library(Seurat)
	library(dplyr)
	library(harmony)
        library(future)
        library(future.apply)
})

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
outfile <- args[2]
nFeature <- 3000

nworker <- min(as.numeric(args[3]),length(availableWorkers()))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 70*1024^3)
set.seed(1129L)

normBySCT <- args[4] ##T or F
normMethod <- ifelse(normBySCT=="T", "SCT", "LogNormalize")
rmTRBV <- args[5]

se <- readRDS(file=infile)
DefaultAssay(se) <- "RNA"
se <- DietSeurat(se, assays="RNA")

#se$patient <- gsub("N","",se$orig.ident)
#se$patient <- gsub("T","",se$patient)

split_seurat <- SplitObject(se, split.by = "orig.ident")	

for(i in 1:length(split_seurat)){
   split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = F)
   if(normBySCT=="T"){
	regVars=c('mitoCountRatio', 'nFeature_RNA', 'nCount_RNA', 'percent_mito', 'percent.mt')
	regVars=intersect(regVars,colnames(se[[]]))
	load(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/resources/cycle.rda")	
	tryCatch({
            split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], s.features = s_genes, g2m.features = g2m_genes)
            split_seurat[[i]]$cc_difference <- split_seurat[[i]]$S.Score - split_seurat[[i]]$G2M.Score
            split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c(regVars, 'cc_difference'))
	}, error=function(cond){
            cat(sprintf("No cell cycle scoring due to %s",cond))
	    split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = regVars)
        })

   }else{
	split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]], selection.method = "vst", nfeatures = nFeature)	
   }
}

#Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = nFeature)

se <- merge(split_seurat[[1]], y=split_seurat[2:length(split_seurat)], merge.data=T)
VariableFeatures(se) <- integ_features
assay2use <- ifelse(normBySCT=="T","SCT","RNA")

if(rmTRBV=="T"){
	hvg <- se[[assay2use]]@var.features
        hvg <- hvg[!grepl("TRBV",hvg)]
        se <- RunPCA(object=se, features=hvg, verbose=F)  ##https://satijalab.org/seurat/v3.2/sctransform_vignette.html
}else{
        se <- RunPCA(object=se, verbose=F)
}
se <- RunHarmony(object=se, assay.use = assay2use, reduction = "pca", dims.use = 1:50, group.by.vars = "orig.ident", plot_convergence = FALSE)
saveRDS(se, file=outfile)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
