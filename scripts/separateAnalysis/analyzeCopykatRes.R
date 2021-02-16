library(Seurat)
library(data.table)
library(ggplot2)
se <- readRDS(file="./data/GEJ_QCed_sctNorm_BatchCCA_clustStab/EpithelialCells_sct_hmy50_k100_res0.05.rds")
pred <- fread(file="./CNV/ACGEJ_copykat/ACGEJ_copykat_ep_prediction.txt",header=T,data.table=F,stringsAsFactors=F)
se <- subset(se, cells=pred$cell.names)
se$ploidy <- pred$copykat.pred[match(colnames(se),pred$cell.names)]
se$sample <- ifelse(grepl("N",se$orig.ident), "normal", "tumor")
table(se$sample,se$ploidy)
p1 <- DimPlot(se, group.by="ploidy", shuffle=T, seed=1129L, label=F, reduction="umap")
p2 <- DimPlot(se, group.by="sample", shuffle=T, seed=1129L, label=F, reduction="umap")
ggsave(plot=p1/p2, filename="./plot/cellident/escc_EpithelialCells_sct_hmy50_k100_res0.05_copykat.pdf, width = 5, height = 8, units="in",dpi=300)



#https://github.com/navinlabcode/copykat

library(Seurat)
library(copykat)
library(data.table)

args  <- commandArgs(trailingOnly=T)
prefix <- args[1] #including the absolute path
ncpu <- as.numeric(args[2])

pred <- fread(file=paste0(prefix,"_prediction.txt",header=T,data.table=F,stringsAsFactors=F)
CNA <- fread(file=paste0(prefix,"_CNA_results.txt",header=T,data.table=F,stringsAsFactors=F)

pred=pred%>%mutate(celltype=ifelse(grepl("-E_",cell.names),"+","-"))
table(pred$copykat.pred,pred$celltype)
pred=pred%>%filter(celltype=="+")%>%select(-celltype)

tumor.cells <- pred$cell.names[which(pred$copykat.pred=="aneuploid")]
tumor.mat <- CNA[, which(colnames(CNA) %in% tumor.cells)]
hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads=ncpu, method = "euclidean"), method = "ward.D2")
hc.umap <- cutree(hcc,2)

rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
cells <- rbind(subpop,subpop)

heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads = ncpu, 
	    method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

  legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')

