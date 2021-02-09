suppressPackageStartupMessages({
        library(Seurat)
        library(ggplot2)
        library(patchwork)
})

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
outprefix <- args[2] #including absolute path


####DotPlot for ACGEJ

#immune (CD45+,PTPRC), epithelial/cancer (EpCAM+,EPCAM), and stromal (CD10+,MME,fibo or CD31+,PECAM1,endo) 
genes_to_check = c("PTPRC","EPCAM",'PECAM1','MME',"CD3G","CD3E", "CD79A")

imm <- c(1,2,3,4,5,6,9,10,11,12,14,15)
epi <- 0
stro <- c(7,8,13)
se@meta.data$celltype <- ifelse(se@meta.data$integrated_snn_res.0.4  %in% imm ,'immune',
    ifelse(se@meta.data$integrated_snn_res.0.4  %in% epi ,'epi','stromal'))
table(se@meta.data$celltype)
se@meta.data$group=paste(se@meta.data$celltype,se@meta.data$integrated_snn_res.0.4)
p <- DotPlot(se, features = genes_to_check, assay="RNA", group.by='group')
ggsave(plot=p, filename="./plot/markerGenes/ACGEJ_generalMarker_dotPlot.pdf", width=8, height=5, units="in", dpi=300)


#Immune Markers, used for differing between primary immune cell types
cells.use <- row.names(se@meta.data)[which(se@meta.data$celltype=='immune')]
se <- subset(se, cells=cells.use)

tc <- c(2,3,5,6,10,11,14)
bc <- c(1,9)
se@meta.data$celltype <- ifelse(se@meta.data$integrated_snn_res.0.4 %in% tc ,'TCells',
    ifelse(se@meta.data$integrated_snn_res.0.4 %in% bc, 'BCells', 
	ifelse(se@meta.data$integrated_snn_res.0.4 == '4', 'Macrophages',
	ifelse(se@meta.data$integrated_snn_res.0.4 == '12', 'Mast',
	ifelse(se@meta.data$integrated_snn_res.0.4 == '15', 'Dendritic', 'Unknown')))))
se@meta.data$group=paste(se@meta.data$celltype,se@meta.data$integrated_snn_res.0.4)
genes_to_check=fread(file="./resources/broad_cell_markers_immune.csv",header=T,data.table=F)
genes_to_check=setdiff(unique(genes_to_check[genes_to_check$cell!="Neutrophils", ]$gene),'CCL3L1')
# genes_to_check=unique(c( 'CD2','CD3D','CD3E','CD3G','MARCO','CSF1R','CD68','GLDN','APOE','CCL3L1',
  # 'TREM2','C1QB','NUPR1','FOLR2','RNASE1','C1QA','CD1E','CD1C','FCER1A','PKIB',
  # 'CYP2S1','NDRG2','CMA1','MS4A2','TPSAB1','TPSB2','IGLL5','MZB1','JCHAIN','DERL3',
  # 'SDC1','MS4A1','BANK1','PAX5','CD79A','PRDM1','XBP1','IRF4','MS4A1','IRF8','ACTB',
  # 'GAPDH','MALAT1','FCGR3B','ALPL','CXCR1','CXCR2','ADGRG3','CMTM2','PROK2','MME','MMP25',
  # 'TNFRSF10C','SLC32A1','SHD','LRRC26','PACSIN1','LILRA4','CLEC4C','DNASE1L3',
  # 'CLEC4C','LRRC26','SCT','LAMP5'))
p <- DotPlot(se, features = genes_to_check, assay="RNA", group.by='group') #+ coord_flip()
p + theme(axis.text.x=element_text(size=0.5,angle=90),axis.text.y=element_text(size=0.5,color="blue"),legend.text=element_text(size=0.5))
ggsave(plot=p, filename="./plot/markerGenes/ACGEJ_immuneMarker_dotPlot.pdf", width=20, height=5, units="in", dpi=300)


####DotPlot for ESCC

genes_to_check = c("PTPRC","EPCAM",'PECAM1','MME',"CD3G","CD3E", "CD79A","KRT8","KRT18")
imm <- c(1,3,9,11,5,10,4,15,14)
epi <- c(0,12)
stro <- c(2,6,7,13,16)
others <- 8
se@meta.data$celltype <- ifelse(se@meta.data$SCT_snn_res.0.2  %in% imm ,'immune',
    ifelse(se@meta.data$SCT_snn_res.0.2  %in% epi ,'epi',
    ifelse(se@meta.data$SCT_snn_res.0.2  %in% stro, 'stromal', 'others')))
table(se@meta.data$celltype)
se@meta.data$group=paste(se@meta.data$celltype,se@meta.data$SCT_snn_res.0.2)
p <- DotPlot(se, features = genes_to_check, assay="RNA", group.by='group')
ggsave(plot=p, filename="./plot/markerGenes/ESCC_generalMarker_dotPlot.pdf", width=10, height=5, units="in", dpi=300)

#General Cell Markers, used for differing between non-immune and immune cell types as well as non immuen epithelial cell types
genes_to_check=c('PTPRC','CD3G','CD3E','CD79A','BLNK','CD68','CSF1R','MARCO','CD207','PMEL','MLANA','PECAM1','CD34','VWF','EPCAM','SFN','KRT19','ACTA2','MCAM','MYLK','MYL9','FAP','THY1','ALB')
se@meta.data$celltype2=paste(se@meta.data$celltype,se@meta.data$SCT_snn_res.0.2)

#Immune Markers, used for differing between primary immune cell types
cells.use <- row.names(se@meta.data)[which(se@meta.data$celltype=='immune')]
se <- subset(se, cells=cells.use)
tc <- c(1,3,9,11)
bc <- c(5,10)
se@meta.data$celltype <- ifelse(se@meta.data$SCT_snn_res.0.2 %in% tc ,'TCells',
    ifelse(se@meta.data$SCT_snn_res.0.2 %in% bc, 'BCells', 
	ifelse(se@meta.data$SCT_snn_res.0.2 == '4', 'Macrophages',
	ifelse(se@meta.data$SCT_snn_res.0.2 == '14', 'Mast',
	ifelse(se@meta.data$SCT_snn_res.0.2 == '15', 'Dendritic', 'Unknown')))))
se@meta.data$group=paste(se@meta.data$celltype,se@meta.data$SCT_snn_res.0.2)
genes_to_check=fread(file="./resources/broad_cell_markers_immune.csv",header=T,data.table=F)
genes_to_check=setdiff(unique(genes_to_check[genes_to_check$cell!="Neutrophils", ]$gene),'CCL3L1')
p <- DotPlot(se, features = genes_to_check, assay="RNA", group.by='group') #+ coord_flip()
p + theme(axis.text.x=element_text(size=0.5,angle=90),axis.text.y=element_text(size=0.5,color="blue"),legend.text=element_text(size=0.5))
ggsave(plot=p, filename="./plot/markerGenes/ESCC_immuneMarker_dotPlot.pdf", width=20, height=5, units="in", dpi=300)




#markergenes=c("MB21D1","TMEM173", "NFKB1", "IFI16", "IFI44", "IFI30", "ITGAM","MAFB", "CSF2", "IL1A","IL1B", "IFNG","TNF", "CCL2","CCL5")
#markergenes=c("CD14","CCL22","FOXP3","CCL17","IL10")
#checkpoints=c("PDCD1","CD274","PDCD1LG2","CTLA4","HAVCR2")
#markergenes=c("TGFB1","CD80","CD86","IDO1","CCL8","CXCL10","FCGR3A")
#markergenes=c("CD3E","CD3D","CD79A")
#markergenes=c("BCL2L1","BCL2","BAX","MCL1","MYC","MYB","BAK1","CASP9","IFNB1")
#markergenes <- c("TLR7","TLR8","TRIM22","CCR5","CD247","HAVCR2","STAT1","FOXP3","LAG3","IDO1","IL10")
#markergenes=c(,,,) #dsRNA sensor genes: MDA5 > IFIH1, RIG-1 > DDX58, LGP2 > DHX58, PKR > EIF2AK2, PACT > PRKRA, TRBP > TARBP2, DDX3 > DDX3X

#markergenes <- list()
#markergenes[[1]] <- c("NLRP1","IFIH1","DDX58","DHX58","EIF2AK2")
#markergenes[[2]] <- c("OAS1","OAS2","OAS3","OASL","ADAR")
#markergenes[[3]] <- c("PRKRA","TARBP2","DHX9","DDX60","DHX15")
#markergenes[[4]] <- c("DHX33","DHX29","DDX3X","DDX17","DDX1","DDX21","DHX36")

se=readRDS(file=infile)
allgenes=rownames(se)
print(paste(intersect(markergenes,allgenes),collapse=","))
markergenes=intersect(markergenes,allgenes)

#for(i in 1:length(markergenes)){
#    jpeg(file=paste0("plot/cellident/ESCC_IM_IFI30corr",i,".jpg"),width=40,height=ceiling(length(markergenes[[i]])/2)*15,unit="cm",res=300)
#    FeaturePlot(se,reduction="glmpca20_umap",features=markergenes[[i]],order=T,min.cutoff='q10',repel=T,ncol=min(length(markergenes[[i]]),2))
#    dev.off()
#}

p <- FeaturePlot(se,reduction="glmpca10_umap",features=markergenes,order=T,min.cutoff='q10',repel=T,ncol=min(length(markergenes),2))
ggsave(plot=p, filename=paste0(outprefix, "_feaPlot.pdf"), width=3, height=ceiling(length(markergenes[[i]])/2)), units="in", dep=300)

p <- DotPlot(se, features = markergenes) + coord_flip()
ggsave(plot=p, filename=paste0(outprefix, "_dotPlot.pdf"), width=5, height=4, units="in", dpi=300)
