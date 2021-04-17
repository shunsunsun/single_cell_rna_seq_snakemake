##Ref: https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
##Ref: https://satijalab.org/seurat/articles/sctransform_vignette.html

options(stringsAsFactors = F)
suppressPackageStartupMessages({
        library(sctransform)
        library(Seurat)
	library(dplyr)
        library(future)
        library(future.apply)
})

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
outfile <- args[2]
nFeature <- 3000
mode <- "cca"

# prevent small datasets from not having enough neighbors (k) to use when filtering anchors
# https://github.com/satijalab/seurat/issues/997#issuecomment-505207415
k_filter <- 200

nworker <- min(as.numeric(args[3]),length(availableWorkers()))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 70*1024^3)
set.seed(1129L)

normBySCT <- args[4] ##T or F
normMethod <- ifelse(normBySCT=="T", "SCT", "LogNormalize")

#if(grepl("TCell", infile) && grepl("ESCC", infile)){
#	refs <- c("GSM4317409_S133A","GSM4317412_S149A","GSM4317413_S150A") #For ESCC_jiaoda_TCells
#}
#if(grepl("TCell", infile) && grepl("GEJ", infile)){
#	refs <- c("LZE8T","LZE7T","LZE6N","LZE8D1") #For ESCC_14TDN_TCells
#}
#if(grepl("BCell", infile) && grepl("ESCC", infile)){
#	refs <- c("LZE11N","LZE22D","LZE3T") #For ESCC_14TDN_BCells
#}
#if(grepl("BCell", infile) && grepl("jiaoda", infile)){
#	refs <- c("GSM4317411_S135A","GSM4317421_S158B") #For ESCC_jiaoda_BCells
#}

se <- readRDS(file=infile)
DefaultAssay(se) <- "RNA"
se <- DietSeurat(se, assays="RNA")

#se$patient <- gsub("N","",se$orig.ident)
#se$patient <- gsub("T","",se$patient)

#https://github.com/satijalab/seurat/issues/3207#issuecomment-651862274
meta <- se@meta.data
tmp<-as.data.frame(meta%>%group_by(orig.ident)%>%summarize(count=n())%>%ungroup())
cells.use <- rownames(meta[meta$orig.ident %in% tmp[tmp$count>k_filter,]$orig.ident,])	
se <- subset(se,cells=cells.use)

split_seurat <- SplitObject(se, split.by = "orig.ident")	

for(i in 1:length(split_seurat)){
   split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
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

if(normBySCT=="T"){
	# Prepare the SCT list object for integration
	split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features)
}

#reference_datasets <- which(names(split_seurat) %in% refs)
reference_datasets <- NULL

# Find best buddies - can take a while to run
if(mode=="cca"){
        integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = normMethod,
                anchor.features = integ_features, reference = reference_datasets, verbose = TRUE, k.filter=k_filter)
}else{
        split_seurat <- future_lapply(X = split_seurat, FUN = RunPCA, verbose=FALSE, features=integ_features)
        integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = normMethod,
                anchor.features = integ_features, reference = reference_datasets, reduction = "rpca", verbose = TRUE, k.filter=k_filter)
}

# Integrate across batches and get a batch-corrected assay of all features (genes)
# https://github.com/satijalab/seurat/issues/2590#issuecomment-622463806
all_features <- Reduce(intersect, lapply(split_seurat, rownames))
se <- IntegrateData(anchorset = integ_anchors, normalization.method = normMethod, features.to.integrate = all_features, verbose=T)
saveRDS(se, file=outfile)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
