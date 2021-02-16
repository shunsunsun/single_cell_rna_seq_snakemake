library(Seurat)
library(copykat)

args  <- commandArgs(trailingOnly=T)
ncpu <- as.numeric(args[1])
dist <- "euclidean"
#infile <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/CNV/ACGEJ_copykat/ACGEJ_copykatInput_clust0478.rds"
#infile <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/CNV/ACGEJ_copykat/ACGEJ_copykatInput_allEpi.rds"
infile <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/CNV/ESCC_copykat/ESCC_copykatInput_allEpi.rds"
#outprefix <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/CNV/ACGEJ_copykat/ACGEJ_copykat_EpithelialCells_"

exp.rawdata <- readRDS(infile)
copykat.res <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.15, distance=dist, n.cores=ncpu)
pred <- data.frame(copykat.res$prediction)
CNA <- data.frame(copykat.res$CNAmat)
#saveRDS(pred, paste0(outprefix, "pred.rds"))
#saveRDS(CNA, paste0(outprefix, "cna.rds"))
