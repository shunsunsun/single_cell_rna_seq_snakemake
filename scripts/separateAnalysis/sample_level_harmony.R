options(stringsAsFactors = F)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(harmony))

args  <- commandArgs(trailingOnly=T)
#celltype <- args[1] #TCells, BCells, Monocytes, ImmuneCells, EpithelialCells, StromalCells
#workdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/Combined/"
#infile <- paste0(workdir,celltype,"_sctNorm_CCA.rds")
#targetAssay <- "integrated"
#outfile <- gsub(".rds","_Hmy.rds",infile)

#infile <- "../../data/ESCC_QCed_sctNorm_BatchHmy_clustStab/EpithelialCells_sctNorm.rds"
infile <- args[1]
#outfile <- "../../data/ESCC_QCed_sctNorm_BatchHmy_clustStab/EpithelialCells_sctNorm_BatchHmy.rds"
outfile <- args[2]
rmTRBV <- args[3] #T or F
targetAssay <- "SCT"


# https://github.com/immunogenomics/harmony/issues/41#issuecomment-633885490
# if samples from different techical platforms, try https://github.com/immunogenomics/harmony/issues/41#issuecomment-642862186

se <- readRDS(file=infile)
se@meta.data$orig.ident=sapply(strsplit(se@meta.data$orig.ident,"-"),"[",1)
if(rmTRBV=="T"){
	hvg <- se@assays$SCT@var.features
	hvg <- hvg[!grepl("TRBV",hvg)]
	se <- RunPCA(object=se, features=hvg, verbose=T)  ##https://satijalab.org/seurat/v3.2/sctransform_vignette.html
}else{
	se <- RunPCA(object=se, verbose=T)
}
se <- RunHarmony(object=se, assay.use = targetAssay, reduction = "pca", dims.use = 1:50, group.by.vars = "orig.ident", plot_convergence = FALSE)
saveRDS(se, file=outfile)
