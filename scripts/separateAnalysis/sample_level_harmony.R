options(stringsAsFactors = F)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(harmony))

args  <- commandArgs(trailingOnly=T)
celltype <- args[1] #TCells, BCells, Monocytes, ImmuneCells, EpithelialCells, StromalCells
workdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/Combined/"
infile <- paste0(workdir,celltype,"_sctNorm_CCA.rds")
targetAssay <- "integrated"
outfile <- gsub(".rds","_Hmy.rds",infile)

# https://github.com/immunogenomics/harmony/issues/41#issuecomment-633885490
# if samples from different techical platforms, try https://github.com/immunogenomics/harmony/issues/41#issuecomment-642862186

se <- readRDS(file=infile)
se@meta.data$orig.ident=sapply(strsplit(se@meta.data$orig.ident,"-"),"[",1)
se <- RunPCA(object=se, verbose=T)  ##https://satijalab.org/seurat/v3.2/sctransform_vignette.html
se <- RunHarmony(object=se, assay.use = targetAssay, reduction = "pca", dims.use = 1:50, group.by.vars = "orig.ident", plot_convergence = FALSE)
saveRDS(se, file=outfile)
