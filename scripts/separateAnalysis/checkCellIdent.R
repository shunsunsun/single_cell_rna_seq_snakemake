library(Seurat)
se <- readRDS(file="./data/GEJ_QCed_sctNorm_BatchCCA_clustStab/Louvain_clust_k100_pc50_res0.4.rds")
load(file="./data/GEJ_QCed_sctNorm_BatchCCA_clustStab/allcells_ident_scibet.rda")
load(file="./data/GEJ_QCed_sctNorm_BatchCCA_clustStab/Louvain_clust_k100_pc50_res0.4_singleR_label.fine.rda")
se$singleR_label <- clust.pred$labels[match(se$integrated_snn_res.0.4, rownames(clust.pred))]
se$scibet_label <- ci[[2]]
table(se$scibet_label, se$singleR_label)
table(se$scibet_label, se$integrated_snn_res.0.4)


