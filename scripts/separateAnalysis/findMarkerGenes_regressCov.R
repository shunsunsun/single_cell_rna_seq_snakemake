#Ref: https://hbctraining.github.io/scRNA-seq/lessons/09_merged_SC_marker_identification.html
#https://github.com/satijalab/seurat/issues/2180

suppressPackageStartupMessages({
	library(Seurat)
	library(future)
	library(tidyverse)
	library(MAST)
})

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
sctransform <- args[2] #T or F
assay_used <- ifelse(sctransform=="T","SCT","RNA") #RNA is recommended; I only find SCT significantly different from RNA results and useful for integrated datasets
affix <- ifelse(sctransform=="T","_sct","")
affix <- paste0(affix, "_regCov")
outfile <- gsub(".rds",paste0(affix,"_allMarkers.rda"),infile)
annotations <- read.csv("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/resources/annotation.csv")

nworker <- min(as.numeric(args[3]),length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3) #60G
randSeed <- 1129L
set.seed(randSeed)

se <- readRDS(file=infile)
clustRes <- colnames(se@meta.data)[grepl("_res\\.", colnames(se@meta.data))]
Idents(se) <- clustRes
#DefaultAssay(se) <- "RNA"

all_markers <- FindAllMarkers(object = se, assay=assay_used, slot="data", only.pos = TRUE, 
	min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST", latent.vars="orig.ident")
save(all_markers, file=outfile)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)


##Analyze the results
##Check analyzeMarkerGenes.R
