library(Seurat)
#se <- readRDS(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/EpithelialCells.rds")
#se <- readRDS(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustStab/Louvain_clust_k100_pc50_res0.4.rds")
#meta <- se@meta.data
#clustCol <- which(grepl("_res\\.",colnames(meta)))
#colnames(meta)[clustCol]="clustRes"
#cells.use <- rownames(meta[meta$clustRes %in% c(0,4,7,8),])
#se <-subset(se, cells=cells.use)

#se <- readRDS(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustStab/EpithelialCells.rds")
se <- readRDS(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/EpithelialCells.rds")
exp.rawdata <- as.matrix(se@assays$RNA@counts)
#saveRDS(exp.rawdata, file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/CNV/ACGEJ_copykat/ACGEJ_copykatInput_allEpi.rds")
saveRDS(exp.rawdata, file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/CNV/ESCC_copykat/ESCC_copykatInput_allEpi.rds")
