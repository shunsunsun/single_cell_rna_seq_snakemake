suppressPackageStartupMessages({
    library(Seurat)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(dplyr)
    library(ggplot2)
    library(tidyverse)
    library(patchwork)
}) 
  
plotdir <- '/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/markerGenes/'
cancertype <- 'ESCC'
infile="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/TCells_rc_glmpca20_k100_res0.15.rds"
se <- readRDS(file=infile)
load(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/TCells_rc_glmpca20_k100_res0.15_allMarkers.rda")
cl=3
sig_dge.celltype <- as.data.frame(all_markers %>% filter(!grepl("^MT-",gene)) %>% filter(!grepl("^RP",gene)) %>% filter(p_val_adj<=0.05, cluster==cl))
sig_dge.celltype <- as.data.frame(all_markers %>% filter(p_val_adj<=0.05, cluster==cl))
sig_dge.celltype <- as.data.frame(all_markers %>% filter(p_val_adj<=0.05, cluster==cl, abs(avg_logFC)>1))

##For T cells
sig_dge.celltype[sig_dge.celltype$gene=="GLTSCR2", ]$gene="NOP53"
sig_dge.celltype[sig_dge.celltype$gene=="GNB2L1", ]$gene="RACK1"
sig_dge.celltype[sig_dge.celltype$gene=="CCDC109B", ]$gene="MCUB"
sig_dge.celltype[sig_dge.celltype$gene=="AIM1", ]$gene="CRYBG1"
sig_dge.celltype[sig_dge.celltype$gene=="FYB", ]$gene="FYB1"
sig_dge.celltype[sig_dge.celltype$gene=="RARRES3", ]$gene="PLAAT4"
sig_dge.celltype[sig_dge.celltype$gene=="TMEM173", ]$gene="STING1"
sig_dge.celltype[sig_dge.celltype$gene=="SEPT6", ]$gene="SEPTIN6"

##For B cells
sig_dge.celltype[sig_dge.celltype$gene=="LRMP", ]$gene="IRAG2"
sig_dge.celltype[sig_dge.celltype$gene=="H2AFZ", ]$gene="H2AZ1"
sig_dge.celltype[sig_dge.celltype$gene=="H2AFV", ]$gene="H2AZ2"
sig_dge.celltype[sig_dge.celltype$gene=="ATP5L", ]$gene="ATP5MG"
sig_dge.celltype[sig_dge.celltype$gene=="ATP5G3", ]$gene="ATP5MC3"
sig_dge.celltype[sig_dge.celltype$gene=="HN1", ]$gene="JPT1"
sig_dge.celltype[sig_dge.celltype$gene=="H2AFY", ]$gene="MACROH2A1"
sig_dge.celltype[sig_dge.celltype$gene=="KIAA0922", ]$gene="TMEM131L"
sig_dge.celltype[sig_dge.celltype$gene=="KIAA0101", ]$gene="PCLAF"
sig_dge.celltype[sig_dge.celltype$gene=="H3F3A", ]$gene="H3-3A"
sig_dge.celltype[sig_dge.celltype$gene=="ATP5J2", ]$gene="ATP5MF"
sig_dge.celltype[sig_dge.celltype$gene=="ATP5C1", ]$gene="ATP5F1C"
sig_dge.celltype[sig_dge.celltype$gene=="TROVE2", ]$gene="RO60"
sig_dge.celltype[sig_dge.celltype$gene=="FAM129A", ]$gene="NIBAN1"
sig_dge.celltype[sig_dge.celltype$gene=="C14orf166", ]$gene="RTRAF"
sig_dge.celltype[sig_dge.celltype$gene=="ATPIF1", ]$gene="ATP5IF1"
sig_dge.celltype[sig_dge.celltype$gene=="ATP5F1", ]$gene="ATP5PB"
sig_dge.celltype[sig_dge.celltype$gene=="SHFM1", ]$gene="SEM1"
sig_dge.celltype[sig_dge.celltype$gene=="FAM49B", ]$gene="CYRIB"
sig_dge.celltype[sig_dge.celltype$gene=="ATP5A1", ]$gene="ATP5F1A"
sig_dge.celltype[sig_dge.celltype$gene=="ATP5J", ]$gene="ATP5PF"
sig_dge.celltype[sig_dge.celltype$gene=="ATP5B", ]$gene="ATP5F1B"
sig_dge.celltype[sig_dge.celltype$gene=="SELT", ]$gene="SELENOT"
sig_dge.celltype[sig_dge.celltype$gene=="C12orf49", ]$gene="SPRING1"
sig_dge.celltype[sig_dge.celltype$gene=="CASC5", ]$gene="KNL1"
sig_dge.celltype[sig_dge.celltype$gene=="SGK223", ]$gene="PRAG1"
sig_dge.celltype[sig_dge.celltype$gene=="FAM195B", ]$gene="MCRIP1"
sig_dge.celltype[sig_dge.celltype$gene=="CCDC109B", ]$gene="MCUB"
sig_dge.celltype[sig_dge.celltype$gene=="SEPT2", ]$gene="SEPTIN2"
sig_dge.celltype[sig_dge.celltype$gene=="TCEB2", ]$gene="ELOB"
sig_dge.celltype[sig_dge.celltype$gene=="H2AFX", ]$gene="H2AX"
sig_dge.celltype[sig_dge.celltype$gene=="ZCCHC6", ]$gene="TUT7"
sig_dge.celltype[sig_dge.celltype$gene=="C11orf31", ]$gene="SELENOH"
sig_dge.celltype[sig_dge.celltype$gene=="DIRC2", ]$gene="SLC49A4"
sig_dge.celltype[sig_dge.celltype$gene=="C14orf2", ]$gene="ATP5MPL"
sig_dge.celltype[sig_dge.celltype$gene=="ATP5E", ]$gene="ATP5F1E"
sig_dge.celltype[sig_dge.celltype$gene=="FAM96B", ]$gene="CIAO2B"
sig_dge.celltype[sig_dge.celltype$gene=="TCEB1", ]$gene="ELOC"
sig_dge.celltype[sig_dge.celltype$gene=="CPSF3L", ]$gene="INTS11"
sig_dge.celltype[sig_dge.celltype$gene=="C7orf73", ]$gene="STMP1"
sig_dge.celltype[sig_dge.celltype$gene=="RTFDC1", ]$gene="RTF2"
sig_dge.celltype[sig_dge.celltype$gene=="GARS", ]$gene="GARS1"
sig_dge.celltype[sig_dge.celltype$gene=="FAM60A", ]$gene="SINHCAF"
sig_dge.celltype[sig_dge.celltype$gene=="C11orf73", ]$gene="HIKESHI"
sig_dge.celltype[sig_dge.celltype$gene=="MINOS1", ]$gene="MICOS10"
sig_dge.celltype[sig_dge.celltype$gene=="MGEA5", ]$gene="OGA"
sig_dge.celltype[sig_dge.celltype$gene=="C20orf24", ]$gene="RAB5IF"
sig_dge.celltype[sig_dge.celltype$gene=="FAM96A", ]$gene="CIAO2A"
sig_dge.celltype[sig_dge.celltype$gene=="SEPT1", ]$gene="SEPTIN1"
sig_dge.celltype[sig_dge.celltype$gene=="SEP15", ]$gene="SELENOF"
sig_dge.celltype[sig_dge.celltype$gene=="SEPT9", ]$gene="SEPTIN9"
sig_dge.celltype[sig_dge.celltype$gene=="DARS", ]$gene="DARS1"
sig_dge.celltype[sig_dge.celltype$gene=="C19orf43", ]$gene="TRIR"
sig_dge.celltype[sig_dge.celltype$gene=="USMG5", ]$gene="ATP5MD"
sig_dge.celltype[sig_dge.celltype$gene=="WHSC1L1", ]$gene="NSD3"
sig_dge.celltype[sig_dge.celltype$gene=="HIST1H4C", ]$gene="H4C3"
sig_dge.celltype[sig_dge.celltype$gene=="HIST1H1E", ]$gene="H1-4"
sig_dge.celltype[sig_dge.celltype$gene=="FAM208B", ]$gene="TASOR2"
sig_dge.celltype[sig_dge.celltype$gene=="C20orf196", ]$gene="SHLD1"
sig_dge.celltype[sig_dge.celltype$gene=="SEPT7", ]$gene="SEPTIN7"
sig_dge.celltype[sig_dge.celltype$gene=="FAM103A1", ]$gene="RAMAC"
sig_dge.celltype[sig_dge.celltype$gene=="MARCH5", ]$gene="MARCHF5"
sig_dge.celltype[sig_dge.celltype$gene=="QARS", ]$gene="QARS1"
sig_dge.celltype[sig_dge.celltype$gene=="RFWD2", ]$gene="COP1"
sig_dge.celltype[sig_dge.celltype$gene=="TMEM261", ]$gene="DMAC1"
sig_dge.celltype[sig_dge.celltype$gene=="MARCH7", ]$gene="MARCHF7"
sig_dge.celltype[sig_dge.celltype$gene=="WDR34", ]$gene="DYNC2I2"
sig_dge.celltype[sig_dge.celltype$gene=="C14orf142", ]$gene="GON7"
sig_dge.celltype[sig_dge.celltype$gene=="TSSC1", ]$gene="EIPR1"
sig_dge.celltype[sig_dge.celltype$gene=="ATP5S", ]$gene="ATP5F1E"
sig_dge.celltype[sig_dge.celltype$gene=="C18orf8", ]$gene="RMC1"
sig_dge.celltype[sig_dge.celltype$gene=="FAM192A", ]$gene="PSME3IP1"
sig_dge.celltype[sig_dge.celltype$gene=="MKL1", ]$gene="MRTFA"
sig_dge.celltype[sig_dge.celltype$gene=="ERBB2IP", ]$gene="ERBIN"
sig_dge.celltype[sig_dge.celltype$gene=="KIAA1033", ]$gene="WASHC4"
sig_dge.celltype[sig_dge.celltype$gene=="EBLN3", ]$gene="EBLN3P"
sig_dge.celltype[sig_dge.celltype$gene=="YARS", ]$gene="YARS1"
sig_dge.celltype[sig_dge.celltype$gene=="ADPRHL2", ]$gene="ADPRS"
sig_dge.celltype[sig_dge.celltype$gene=="MATR3.1", ]$gene="MATR3"
sig_dge.celltype[sig_dge.celltype$gene=="ASNA1", ]$gene="GET3"
sig_dge.celltype[sig_dge.celltype$gene=="STRA13", ]$gene="BHLHE40"
sig_dge.celltype[sig_dge.celltype$gene=="FOPNL", ]$gene="CEP20"
sig_dge.celltype[sig_dge.celltype$gene=="NUPL2", ]$gene="NUP42"
sig_dge.celltype[sig_dge.celltype$gene=="ADRBK1", ]$gene="GRK2"
sig_dge.celltype[sig_dge.celltype$gene=="C17orf62", ]$gene="CYBC1"
sig_dge.celltype[sig_dge.celltype$gene=="C3orf17", ]$gene="NEPRO"
sig_dge.celltype[sig_dge.celltype$gene=="KIAA0430", ]$gene="MARF1"
sig_dge.celltype[sig_dge.celltype$gene=="KIAA0141", ]$gene="DELE1"
sig_dge.celltype[sig_dge.celltype$gene=="FAM65B", ]$gene="RIPOR2"
sig_dge.celltype[sig_dge.celltype$gene=="FAM129C", ]$gene="NIBAN3"
sig_dge.celltype[sig_dge.celltype$gene=="C6orf48", ]$gene="SNHG32"
sig_dge.celltype[sig_dge.celltype$gene=="IARS", ]$gene="IARS1"
sig_dge.celltype[sig_dge.celltype$gene=="KARS", ]$gene="KARS1"
sig_dge.celltype[sig_dge.celltype$gene=="HARS", ]$gene="HARS1"
sig_dge.celltype[sig_dge.celltype$gene=="C21orf59", ]$gene="CFAP298"
sig_dge.celltype[sig_dge.celltype$gene=="COL4A3BP", ]$gene="CERT1"
sig_dge.celltype[sig_dge.celltype$gene=="C17orf89", ]$gene="NDUFAF8"
sig_dge.celltype[sig_dge.celltype$gene=="ATP5D", ]$gene="ATP5F1D"
sig_dge.celltype[sig_dge.celltype$gene=="VIMP", ]$gene="SELENOS"
sig_dge.celltype[sig_dge.celltype$gene=="LINC00493", ]$gene="SMIM26"
sig_dge.celltype[sig_dge.celltype$gene=="PAPD4", ]$gene="TENT2"
sig_dge.celltype[sig_dge.celltype$gene=="TRAPPC2P1", ]$gene="TRAPPC2B"
sig_dge.celltype[sig_dge.celltype$gene=="MARCH6", ]$gene="MARCHF6"
sig_dge.celltype[sig_dge.celltype$gene=="LINC01420", ]$gene="NBDY"
sig_dge.celltype[sig_dge.celltype$gene=="TSTA3", ]$gene="GFUS"
sig_dge.celltype[sig_dge.celltype$gene=="MARCH1", ]$gene="MARCHF1"
sig_dge.celltype[sig_dge.celltype$gene=="SSSCA1", ]$gene="ZNRD2"
sig_dge.celltype[sig_dge.celltype$gene=="MMP24-AS1", ]$gene="MMP24OS"
sig_dge.celltype[sig_dge.celltype$gene=="HMHA1", ]$gene="ARHGAP45"
sig_dge.celltype[sig_dge.celltype$gene=="ATP5G1", ]$gene="ATP5MC1"
sig_dge.celltype[sig_dge.celltype$gene=="ADSS", ]$gene="ADSS2"
sig_dge.celltype[sig_dge.celltype$gene=="ZCCHC11", ]$gene="TUT4"
sig_dge.celltype[sig_dge.celltype$gene=="TMEM55B", ]$gene="PIP4P1"
sig_dge.celltype[sig_dge.celltype$gene=="EPRS", ]$gene="EPRS1"

#mainly for myeloids
sig_dge.celltype[sig_dge.celltype$gene=="LINC00936", ]$gene="ATP2B1-AS1"
sig_dge.celltype[sig_dge.celltype$gene=="FAM26F", ]$gene="CALHM6"
sig_dge.celltype[sig_dge.celltype$gene=="LINC01272", ]$gene="PELATON"
#sig_dge.celltype[sig_dge.celltype$gene=="FYB1", ]$gene="FYB1"
sig_dge.celltype[sig_dge.celltype$gene=="CECR1", ]$gene="ADA2"
#sig_dge.celltype[sig_dge.celltype$gene=="MACROH2A1", ]$gene="MACROH2A1"
sig_dge.celltype[sig_dge.celltype$gene=="SQRDL", ]$gene="SQOR"
#sig_dge.celltype[sig_dge.celltype$gene=="MCUB", ]$gene="MCUB"
#sig_dge.celltype[sig_dge.celltype$gene=="CYRIB", ]$gene="CYRIB"
#sig_dge.celltype[sig_dge.celltype$gene=="MT-ND4", ]$gene="MT-ND4"
sig_dge.celltype[sig_dge.celltype$gene=="MARCHF1", ]$gene="MARCHF1"
sig_dge.celltype[sig_dge.celltype$gene=="H3-3A", ]$gene="H3-3A"
sig_dge.celltype[sig_dge.celltype$gene=="C10orf54", ]$gene="VSIR"
sig_dge.celltype[sig_dge.celltype$gene=="LINC00152", ]$gene="CYTOR"
#sig_dge.celltype[sig_dge.celltype$gene=="ATP5F1E", ]$gene="ATP5F1E"
#sig_dge.celltype[sig_dge.celltype$gene=="ELOC", ]$gene="ELOC"
sig_dge.celltype[sig_dge.celltype$gene=="WARS", ]$gene="WARS1"

#sig_dge.celltype[sig_dge.celltype$gene=="", ]$gene=""


########
## ORA

#GO
ego_ALL <- enrichGO(gene         = unique(sig_dge.celltype$gene),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_ALL <- data.frame(ego_ALL)
as.data.frame(ego_ALL)%>%arrange(p.adjust)%>%select(Description,GeneRatio)
saveRDS(ego_ALL,file=gsub(".rds",paste0('_C',cl,'marker_GO.rds'),infile))
saveRDS(ego_ALL,file=gsub(".rds",paste0('_C',cl,'marker_rmMTRP_GO.rds'),infile))
       
ego_CC <- enrichGO(gene          = unique(sig_dge.celltype$gene),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.01)
ego_MF <- enrichGO(gene          = unique(sig_dge.celltype$gene),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.01)
ego_BP <- enrichGO(gene          = unique(sig_dge.celltype$gene),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.01)           
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("GO Biological process")
p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("GO Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("GO Molecular function")
plotc <- p_BP/p_CC/p_MF + plot_layout(heights = c(10, 1, 10))
plotc <- p_BP/p_MF + plot_layout(heights = c(20,1))
#plot_layout(widths = c(2, 1), heights = unit(c(5, 1), c('cm', 'null')))
ggsave(paste0(plotdir,cancertype,'_', gsub(".rds", "", basename(infile)),'_C',cl,'marker_GO.pdf'), plotc, width = 12,height = 10)

##KEGG
genes <- bitr(unique(sig_dge.celltype$gene), fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Hs.eg.db')
setdiff(unique(sig_dge.celltype$gene),genes$SYMBOL)
genes <- pull(genes,ENTREZID)               
ekegg <- enrichKEGG(gene = genes, organism = 'hsa', qvalueCutoff=0.05)
p1 <- barplot(ekegg, showCategory=20)
p2 <- dotplot(ekegg, showCategory=20)
plotc = p1/p2
ggsave(paste0(plotdir,cancertype,'_', gsub(".rds", "", basename(infile)),'_C',cl,'marker_KEGG.pdf'), plot = plotc, width = 12, height = 10)

######
# GSEA

#https://github.com/YuLab-SMU/DOSE/wiki/how-to-prepare-your-own-geneList

ranks <- sig_dge.celltype[,c("avg_logFC","gene")]
genes <- bitr(unique(ranks$gene), fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Hs.eg.db')
ranks <- merge(ranks,genes,by.x="gene",by.y="SYMBOL")
ranked.genes <- ranks$avg_logFC
names(ranked.genes) <- ranks$ENTREZID
ranked.genes = sort(ranked.genes, decreasing = TRUE)
 

##GO
gsego_bp <- gseGO(geneList = ranked.genes, OrgDb = org.Hs.eg.db, ont = "BP", nPerm = 5000, pvalueCutoff = 0.1, verbose = FALSE)
gsego_mf <- gseGO(geneList = ranked.genes, OrgDb = org.Hs.eg.db, ont = "MF", nPerm = 5000, pvalueCutoff = 0.1, verbose = FALSE)
gsego_cc <- gseGO(geneList = ranked.genes, OrgDb = org.Hs.eg.db, ont = "CC", nPerm = 5000, pvalueCutoff = 0.1, verbose = FALSE)

##KEGG
gse_kk <- gseKEGG(geneList = ranked.genes, organism = 'hsa', nPerm = 5000, pvalueCutoff = 0.1, verbose = FALSE)


##MsigDB
gmtfile="/home/gaog_pkuhpc/users/liny/geneSets_snps/data/MSigDBgeneset/c5.all.v7.2.entrez.gmt"
term2gene <- read.gmt(gmtfile)
egmt <- enricher(names(ranked.genes), TERM2GENE=term2gene)
head(egmt)
egmt2 <- GSEA(ranked.genes, TERM2GENE=term2gene, verbose=FALSE)
head(egmt2)

##Cell markers
cell_markers <- vroom::vroom('/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/resources/Human_cell_markers.txt') %>%
   tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
   dplyr::select(cellMarker, geneID) %>%
   dplyr::mutate(geneID = strsplit(geneID, ', '))
   
y <- enricher(names(ranked.genes), TERM2GENE=cell_markers, minGSSize=1)
as.data.frame(y)  


  
m_df<- msigdbr(species = "Homo sapiens", category = "C7") #immunological gene sets
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
cluster0.genes<- pbmc.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(p_adj, desc(logFC)) %>%
  dplyr::select(feature, auc)

ranks<- deframe(cluster0.genes)

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) %>%
  head()
  
ggplot(fgseaResTidy %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal() ####以7.5进行绘图填色
  
#比较cluster0和cluster1的差异表达基因
dge.cluster <- FindMarkers(se,ident.1 = 0,ident.2 = 1)
sig_dge.cluster <- subset(dge.cluster, p_val_adj<0.01&abs(avg_logFC)>1)


#比较拟时State1和State3的差异表达基因
p_data <- subset(pData(mycds),select='State')
scRNAsub <- subset(scRNA, cells=row.names(p_data))
scRNAsub <- AddMetaData(scRNAsub,p_data,col.name = 'State')
dge.State <- FindMarkers(scRNAsub, ident.1 = 1, ident.2 = 3, group.by = 'State')
sig_dge.State <- subset(dge.State, p_val_adj<0.01&abs(avg_logFC)>1) 
