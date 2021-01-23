#https://github.com/satijalab/seurat/issues/1883

#Example: sbatch -p cn_icg -A gaog_g1 --qos=gaogcnicg -N 1 -n 1 -c 10 ./runRscript.sh subset_recluster.R \
	# ../../data/GEJ_QCed_sctNorm_BatchCCA_clustStab/Louvain_clust_k100_pc50_res0.4.rds \
	# 2,3,5,6,10,11,14 \
	# ../../data/GEJ_QCed_sctNorm_BatchCCA_clustStab/T_cells.rds


suppressPackageStartupMessages(library(Seurat))

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
tgtClust <- as.character(args[2])
tgtClust <- as.numeric(unname(unlist(strsplit(tgtClust,','))))
outfile <- args[3]

se <- readRDS(file=infile)
meta <- se@meta.data
clustCol <- which(grepl("_res\\.",colnames(meta)))
colnames(meta)[clustCol]="clustRes"
cells.use <- rownames(meta[meta$clustRes %in% tgtClust, ])
se <-subset(se, cells=cells.use)
DefaultAssay(se) <- "RNA"
se <- DietSeurat(se,assays="RNA")
col2rm <- colnames(se@meta.data)[grepl("_SCT|_res\\.",colnames(se@meta.data))]
for(c in col2rm){
	se[[c]] <- NULL
}
saveRDS(se, file=outfile)

