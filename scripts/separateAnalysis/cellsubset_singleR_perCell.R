suppressPackageStartupMessages({
	library(Seurat)
	library(SingleR)
	library(celldex)
	library(BiocParallel)
})

#ref <- "hpca_be_dice_nh_mi"
ref <- "mi"
label <- "label.fine" #label.main or label.fine
refdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/resources/SingleRdata"
samplefile <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/config/sample_cancertype.txt"
randseed=1129L

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
nworkers=min(as.numeric(args[2]),parallel::detectCores())
bpParam=BiocParallel::MulticoreParam(nworkers)

object <- args[3] #se (seurat) or sce
assay <- args[4] #SCT or RNA

if(grepl("Combined", infile)){
	type="combined"
}else if(grepl("ESCC",infile)){
	type="escc"
}else{
	type="gej"
}

outfile <- gsub(".rds","_singleR_perCell_mi.rda",infile)

if(!file.exists(outfile)){
	##--------Get reference---------

	ref_list <- list()
	label_list <- list()
	if(grepl("hpca",ref)){
		ref_list[["hpca"]] <- readRDS(file=paste0(refdir,'/', 'HumanPrimaryCellAtlasData.rds'))
		label_list[["hpca"]] <- ref_list[["hpca"]][[label]]
	}
	if(grepl("be",ref)){
		ref_list[["be"]] <- readRDS(file=paste0(refdir,'/', 'BlueprintEncodeData.rds'))
		label_list[["be"]] <- ref_list[["be"]][[label]]
	}
	if(grepl("dice",ref)){
		ref_list[["dice"]] <- readRDS(file=paste0(refdir,'/', 'DatabaseImmuneCellExpressionData.rds'))
		label_list[["dice"]] <- ref_list[["dice"]][[label]]
	}
	if(grepl("nh",ref)){
		ref_list[["nh"]] <- readRDS(file=paste0(refdir,'/', 'NovershternHematopoieticData.rds'))
		label_list[["nh"]] <- ref_list[["nh"]][[label]]
	}
	if(grepl("mi",ref)){
		ref_list[["mi"]] <- readRDS(file=paste0(refdir,'/', 'MonacoImmuneData.rds'))
		label_list[["mi"]] <- ref_list[["mi"]][[label]]
	}

	##-------Get data------------
	se <- readRDS(file=infile)
	if(object=="sce"){
		se <- as.Seurat(se)
	}
	test=GetAssayData(object = se[[toupper(assay)]], slot = "data")
	##------SingleR identification--------	
	cell.pred <- SingleR(test = test, ref = ref_list, labels = label_list, BPPARAM=bpParam)
	save(cell.pred, file=outfile)
}
