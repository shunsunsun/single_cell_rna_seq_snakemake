suppressPackageStartupMessages({
         library(Seurat)
         library(ggplot2)
         library(patchwork)
})

se <- readRDS(file="data/GEJ_QCed_sctNorm_BatchCCA_clustStab/Louvain_clust_k100_pc50_res0.4.rds")
meta <- se@meta.data
clustResCol <- colnames(meta)[grepl("_res\\.", colnames(meta))]
randseed=1129L
p1 = DimPlot(se, group.by=clustResCol, label=T, reduction="umap") + NoLegend()
p2 = DimPlot(se, group.by="orig.ident", shuffle=T, seed=randseed, label=F, reduction="umap") + NoLegend()
se[["facs"]] <- ifelse(grepl("-E", meta$ident), "CD45-", "CD45+")
p3 = DimPlot(se, group.by="facs", label=F, reduction="umap")
DefaultAssay(se) <- "SCT"
p4 = FeaturePlot(se,reduction="umap",features='PTPRC',order=T,min.cutoff='q10',repel=T, cols = c("green", "blue"))
p = (p1 + p3) / (p2 + p4)
p = p + plot_annotation(tag_levels = 'A')
ggsave(filename="./plot/oralReport_0315/ACGEJ_all_clust.pdf", plot=p, width = 10, height = 8, units="in",dpi=300)
