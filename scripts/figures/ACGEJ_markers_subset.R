##T cells

se <- readRDS(file="data/GEJ_QCed_sctNorm_BatchCCA_clustStab/TCells_sct_hmy20_k100_res0.05.rds")
randseed=1129L
meta <- se@meta.data
clustResCol <- colnames(meta)[grepl("_res\\.", colnames(meta))]
p1 = DimPlot(se, group.by=clustResCol, label=T, reduction="umap") + NoLegend()
p2 = DimPlot(se, group.by="orig.ident", shuffle=T, seed=randseed, label=F, reduction="umap") + NoLegend()

genes_to_check0 <- c("CD4","CD8A","IL2RA","IL7R","FOXP3","TNF","IFNG","KLRK1","NCR1")
p3 <- DotPlot(se, features = genes_to_check0, assay="SCT", group.by=clustResCol) + coord_flip()
p <- (p1+p2)/p3
ggsave(plot=p, filename="./plot/markerGenes/ACGEJ_Tcell.pdf", width=10, height=8, units="in", dpi=300)

#f1 <- fread(file="./resources/T_cell_subset_markers_derLeun.csv",header=T,sep=",",data.table=F)
#f2 <- fread(file="./resources/Tcell_subtypes_literature.csv",header=T,sep=",",data.table=F)
#genes_to_check=unique(union(f1$gene, f2$gene))
#genes_to_check=unique(f1$gene)
#p <- DotPlot(se, features = genes_to_check0, assay="SCT", group.by=clustResCol) + coord_flip()
#ggsave(plot=p, filename="./plot/markerGenes/ACGEJ_Tcell.pdf", width=10, height=8, units="in", dpi=300)
