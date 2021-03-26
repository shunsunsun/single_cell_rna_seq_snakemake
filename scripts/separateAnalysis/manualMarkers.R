gej <- readRDS(file="data/GEJ_QCed_sctNorm_BatchCCA_clustStab/TCells_rc_glmpca10_k150_res0.15_withSCT.rds")
gej$ident <- paste0("Clust",gej$RNA_snn_res.0.15)
levels(gej$ident) <- paste0("Clust",c("0","1","2"))
Idents(gej) <- paste0("Clust",gej$RNA_snn_res.0.15)
levels(gej) <- paste0("Clust",c("0","1","2"))
escc <- readRDS(file="data/ESCC_QCed_sctNorm_BatchHmy_clustStab/TCells_rc_glmpca20_k100_res0.15_withSCT.rds")
escc$ident <- paste0("Clust", escc$RNA_snn_res.0.15)
levels(escc$ident) <- paste0("Clust", c("0","1","2","3"))
Idents(escc)<- paste0("Clust", escc$RNA_snn_res.0.15)
levels(escc) <- paste0("Clust", c("0","1","2","3"))

tcell_id <- c('CD3D','CD3E','CD4','CD8A','CD8B','ZNF683','NCAM1') #NCAM1, CD160 are nk cell markers
#Traditionally, CD3 is never expressed on B or NK cells; ZNF683 -- NKT cells; CD27 -- memory T; CD247 -- encode CD3zeta; CCL5 -- produced by CD8+ T
cytotoxic <- c('GZMB','GNLY','NKG7','FASLG','GZMH','PRF1','GZMA','GZMK','CD40LG')
costimulatory <- c('TNFRSF4','TNFRSF18','TNFRSF1B','TNFRSF9','ICOS','CD27')
coinhibitory <- c('LAG3','CTLA4','HAVCR2','PDCD1','TIGIT','ENTPD1')
survival <- c('CD28','IL32','TMEM173')
location <- c('CCR6','CCR8','ICAM2','ITGA1','ITGAE','ADGRG1')
dys_exhau <- c('NR4A1','NR4A2','NR4A3','RP11-291B21.2')
tf <- c('TCF7','LEF1','ANXA1','RGS2','FOXP3','BATF','RBPJ','TBX21','TOX','EOMES','NFATC2','IRF4','KLF2','SOX4','IKZF2','MIR4435-2HG')
cytokine <- c('CCR7','IL7R','IL2RA','IL1R2','CCL3','CXCL13','CCL4','CCL4L2','CCL5','IFNG','CD55')
migration <- c('SELL','S1PR1')

genesOfinterest <- unique(c(tcell_id,cytotoxic,costimulatory,coinhibitory,survival,location,dys_exhau,tf,cytokine,migration))
p1=DoHeatmap(gej, features = genesOfinterest) + NoLegend() + scale_fill_gradientn(colors=c("blue", "white", "red")) + theme(text = element_text(size = 20))
p2=DoHeatmap(escc, features = genesOfinterest) + NoLegend() + scale_fill_gradientn(colors=c("blue", "white", "red")) + theme(text = element_text(size = 20))
p=p1/p2
ggsave(plot=p,file="./plot/markerGenes/Tcell_heatmap_tf.pdf",width=8, height=10, units="in", dpi=300)

#cytotoxic <- c('CST7','GZMA','GZMB','NKG7','IFNG','PRF1')
#treg <- c('FOXP3','IL2RA')
#naiveT <- c('CCR7','IL7R','LEF1','SELL','TCF7','S1PR1','CD28') ##LEF1 + TCF7 may repress the development of CD4+ lineage
#memoryT <- c('CD27','IL7R','CD44','SELL','CXCR3') ##CD62L -- SELL, CD127 -- IL7R
#reversableExhaust <- c('TBX21','EOMES','TOX') #TBX21 -- T-bet, no EOMES, 
#nkmarkers <- c('NCAM1','CD160')

gej <- readRDS(file="data/GEJ_QCed_sctNorm_BatchCCA_clustStab/BCells_rc_glmpca10_k100_res0.05_withSCT.rds")
gej$ident <- paste0("Clust",gej$RNA_snn_res.0.05)
levels(gej$ident) <- paste0("Clust",c("0","1"))
Idents(gej) <- paste0("Clust",gej$RNA_snn_res.0.05)
levels(gej) <- paste0("Clust",c("0","1"))
escc <- readRDS(file="data/ESCC_QCed_sctNorm_BatchHmy_clustStab/TCells_rc_glmpca10_k150_res0.05_withSCT.rds")
escc$ident <- paste0("Clust", escc$RNA_snn_res.0.05)
levels(escc$ident) <- paste0("Clust", c("0","1","2"))
Idents(escc)<- paste0("Clust", escc$RNA_snn_res.0.05)
levels(escc) <- paste0("Clust", c("0","1","2"))


bcell_id <- c('IGLL5','MZB1','JCHAIN','DERL3','SDC1','BANK1','PAX5','CD79A','CD79B','CD19','PRDM1','XBP1','IRF4','MS4A1','IRF8','CD27','TCL1A')
#MS4A1, IRF8, CD27 -- B mem; PRDM1, XBP1, IRF4 -- Plasma B; naive B -- TCL1A
chemokine <- c('CCR7','CXCR4')
tf <- c('AICDA','BCL6','RGS13','POU2AF1')
bcell_gc <- c('MEF2B','BACH2','MEF2C','EBF1','CD40')
motion <- c('ACTG1','MYO1E','MARCKSL1')
vdj_recomb <- c('HMGB1','HMGB2')
transcription <- c('STMN1','HMGN1','HMGN2')

genesOfinterest <- unique(c(bcell_id,chemokine,tf,bcell_gc,motion,vdj_recomb,transcription))
p1=DoHeatmap(gej, features = genesOfinterest) + NoLegend() + scale_fill_gradientn(colors=c("blue", "white", "red")) + theme(text = element_text(size = 20))
p2=DoHeatmap(escc, features = genesOfinterest) + NoLegend() + scale_fill_gradientn(colors=c("blue", "white", "red")) + theme(text = element_text(size = 20))
p=p1/p2
ggsave(plot=p,file="./plot/markerGenes/Tcell_heatmap_tf.pdf",width=8, height=10, units="in", dpi=300)
