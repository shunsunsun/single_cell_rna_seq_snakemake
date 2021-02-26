imm <- c(1,2,3,4,5,6,9,10,11,12,14,15)
epi <- 0
stro <- c(7,8,13)
se@meta.data$celltype <- ifelse(se@meta.data$integrated_snn_res.0.4  %in% imm ,'immune',
    ifelse(se@meta.data$integrated_snn_res.0.4  %in% epi ,'epi','stromal'))
se@meta.data$group=paste0(se@meta.data$celltype,se@meta.data$integrated_snn_res.0.4)
p <- DotPlot(se, features = genes_to_check, assay="RNA", group.by='group')
ggsave(plot=p, filename="./plot/markerGenes/ACGEJ_generalMarker_dotPlot.pdf", width=8, height=5, units="in", dpi=300)

cells.use <- row.names(se@meta.data)[which(se@meta.data$celltype=='immune')]
se <- subset(se, cells=cells.use)

tc <- c(2,3,5,6,10,11,14)
bc <- c(1,9)
se@meta.data$celltype <- ifelse(se@meta.data$integrated_snn_res.0.4 %in% tc ,'T',
    ifelse(se@meta.data$integrated_snn_res.0.4 %in% bc, 'B', 
	ifelse(se@meta.data$integrated_snn_res.0.4 == '4', 'Macro',
	ifelse(se@meta.data$integrated_snn_res.0.4 == '12', 'Mast',
	ifelse(se@meta.data$integrated_snn_res.0.4 == '15', 'DC', 'Unknown')))))
se@meta.data$group=paste0(se@meta.data$celltype,se@meta.data$integrated_snn_res.0.4)
genes_to_check=fread(file="./resources/broad_cell_markers_immune.csv",header=T,data.table=F)
genes_to_check=setdiff(unique(genes_to_check[genes_to_check$cell!="Neutrophils", ]$gene),'CCL3L1')
p <- DotPlot(se, features = genes_to_check, assay="SCT", group.by='group') + coord_flip()
ggsave(plot=p, filename="./plot/markerGenes/ACGEJ_immuneMarker_dotPlot.pdf", width=10, height=20, units="in", dpi=300)
