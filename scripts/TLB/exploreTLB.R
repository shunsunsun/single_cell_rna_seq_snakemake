suppressPackageStartupMessages({
	library(Seurat)
    library(ggplot2)
    library(patchwork)
})


########
## Process VDJ results

indir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data_vdj"
target_files <- list.files(path=indir,pattern="filtered_contig_annotations.csv",recursive=T)

#It may be worthwhile to remove/ignore columns from your VDJ data that are FALSE or None in the "productive" column, as you'll end up with a lot of extra noise otherwise. Actually, it might be best to just get the clonotype id for each cell, then get the gene information from the consensus_annotations.csv file.
tcr.list <- lapply(target_files, function(f){
        sampleID <- gsub("-VDJ","",dirname(f))
        tcr <- read.csv(file=paste0(indir,"/",f))
        # Remove the -1 at the end of each barcode.
        # Subsets so only the first line of each barcode is kept,
        # as each entry for given barcode will have same clonotype.
        tcr$barcode <- gsub("^", paste0(sampleID,"_"), gsub("-1", "", tcr$barcode))
        tcr <- tcr[!duplicated(tcr$barcode), ]
        # Only keep the barcode and clonotype columns. 
        # We'll get additional clonotype info from the clonotype table.
        tcr <- tcr[,c("barcode", "raw_clonotype_id")]
        colnames(tcr)=gsub("raw_clonotype_id","clonotype_id",colnames(tcr))
        
        # Clonotype-centric info.
        clono <- read.csv(file=paste0(indir,"/",gsub("filtered_contig_annotations","clonotypes", f)))
        
        # Slap the AA sequences onto our original table by clonotype_id.
        tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])

        # Reorder so barcodes are first column and set them as rownames.
        tcr <- tcr[, c(2,1,3)]
        rownames(tcr) <- tcr[,1]
        tcr[,1] <- NULL
        return(tcr)
})

tcr_merged <- do.call(rbind, tcr.list)
gej_samples <- scan(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/config/samples_GEJ.txt",what="character")
gej_samples <- unique(sapply(gej_samples,function(s){
    unname(unlist(strsplit(s,'-')))[1]
}))
tcr_merged$sampleID <- sapply(rownames(tcr_merged), function(r){
    unname(unlist(strsplit(r,'_')))[1]
})
saveRDS(tcr_merged[tcr_merged$sampleID %in% gej_samples, ], file="data/GEJ_TCR_fltContig_clonotype.rds")
saveRDS(tcr_merged[!(tcr_merged$sampleID %in% gej_samples), ], file="data/ESCC_TCR_fltContig_clonotype.rds")

######
## Find TLB

se <- readRDS(file="data/GEJ_IM_pipeCompflted_clustered.rds")
#se <- readRDS(file="data/ESCC_IM_pipeCompflted_clustered.rds")
old_cell_barcodes <- colnames(se)
new_cell_barcodes <- gsub("-I","",gsub("-1","",old_cell_barcodes))
se <- RenameCells(se, new.names=new_cell_barcodes)

allgenes=rownames(se)
markergenes = c('CD79A','CD19','SDC1','IGHG1','IGHA1','JCHAIN')
print(paste(intersect(markergenes,allgenes),collapse=","))
markergenes=intersect(markergenes,allgenes)
jpeg(file="plot/cellident/ESCC_IM_findTLB.jpg",width=40,height=ceiling(length(markergenes)/2)*15,unit="cm",res=300)
FeaturePlot(se,reduction="glmpca20_umap",features=markergenes,order=T,min.cutoff='q10',repel=T,ncol=min(length(markergenes),2))
dev.off()

bcells_in_GEJ_IM <- c(4,7)
bcells_in_GEJ_EP <- 16
bcells_in_ESCC_IM <- c(5,6,8)
bcells_in_ESCC_EP <- c(1,19)

Idents(se) <- 'seurat_clusters'
b_cell <- subset(se, idents = bcells_in_GEJ_IM)
#b_cell <- subset(se, idents = bcells_in_ESCC_IM)
b_cell_new <- CreateSeuratObject(b_cell@assays$RNA@counts ,min.cells = 3, min.genes = 200, project = "gej_b_cell", meta.data = b_cell@meta.data)
#b_cell_new <- CreateSeuratObject(b_cell@assays$RNA@counts ,min.cells = 3, min.genes = 200, project = "escc_b_cell", meta.data = b_cell@meta.data)

b_cell_new <- NormalizeData(object = b_cell_new)
b_cell_new <- FindVariableFeatures(object = b_cell_new)
all.genes <- rownames(b_cell_new)
b_cell_new <- ScaleData(object = b_cell_new, features = all.genes)
b_cell_new <- RunPCA(object = b_cell_new)
b_cell_new <- FindNeighbors(object = b_cell_new)
b_cell_new <- FindClusters(object = b_cell_new)
b_cell_new <- RunUMAP(object = b_cell_new,dims = 1:30)
saveRDS(b_cell_new,file="data/GEJ_IM_TLB.rds")
saveRDS(b_cell_new,file="data/ESCC_IM_TLB.rds")

jpeg(file="plot/cellident/GEJ_Bcells_findTLB.jpg",width=60,height=15,unit="cm",res=300)
p1=DimPlot(object = b_cell_new, reduction = "umap", label = T, cols='alphabet') + NoLegend()
p2=FeaturePlot(b_cell_new, features = c('CD3E', 'CD3D'))
p1+p2+plot_layout(widths = c(20, 40))
dev.off()


#######
# Check whether TLB expressed TCR

escc_tcr <- readRDS(file="data/ESCC_TCR_fltContig_clonotype.rds")
escc_b <- readRDS(file="data/ESCC_IM_TLB.rds")
escc_b_meta <- escc_b@meta.data
escc_tlb=escc_b_meta[escc_b_meta$seurat_clusters %in% c(12, 14),] #605 cells
hit <- escc_tlb[rownames(escc_tlb) %in% rownames(escc_tcr), ] #491 cells
escc_tlb$cellID <- rownames(escc_tlb)
escc_tcr$cellID <- rownames(escc_tcr)
escc_tlb_tcr <- merge(escc_tcr, escc_tlb)

gej_tcr <- readRDS(file="data/GEJ_TCR_fltContig_clonotype.rds")
gej_b <- readRDS(file="data/GEJ_IM_TLB.rds")
gej_b_meta <- gej_b@meta.data
gej_tlb=gej_b_meta[gej_b_meta$seurat_clusters==10,] #225 cells
hit <- gej_tlb[rownames(gej_tlb) %in% rownames(gej_tcr), ] #177 cells
gej_tlb$cellID <- rownames(gej_tlb)
gej_tcr$cellID <- rownames(gej_tcr)
gej_tlb_tcr <- merge(gej_tcr, gej_tlb)


escc_tlb <- escc_tlb[,c(1:4,23,29)]
colnames(escc_tlb)[1]="PatientID"
escc_tlb$TCR <- ifelse(escc_tlb$cellID %in% unique(escc_tcr$cellID), "Y", "N")
escc_tlb$PatientID=gsub("-I","",escc_tlb$PatientID)
rownames(escc_tlb)=NULL
saveRDS(escc_tlb,file="data/ESCC_IM_TLB_meta.rds")

gej_tlb <- gej_tlb[,c(1:4,23,29)]
colnames(gej_tlb)[1]="PatientID"
gej_tlb$TCR <- ifelse(gej_tlb$cellID %in% unique(gej_tcr$cellID), "Y", "N")
gej_tlb$PatientID=gsub("-I","",gej_tlb$PatientID)
rownames(gej_tlb)=NULL
saveRDS(gej_tlb,file="data/GEJ_IM_TLB_meta.rds")


##########
## Association studies

escc_se <- readRDS(file="data/ESCC_IM_TLB.rds")
escc=escc_se@meta.data
colnames(escc)[1]="patientID"
escc$patientID=gsub("-I","",escc$patientID)
escc$TLB=ifelse(escc$seurat_clusters %in% c(12,14),"Y","N")
count <- table(escc$patientID,escc$TLB)
escc <- data.frame(patientID=rownames(count),nTLB=unname(count[,2]), nOtherB=unname(count[,1]))
escc$TLB_B_ratio <- escc$nTLB/escc$nOtherB
escc_adj <- escc[grepl("N", escc$patientID), ]
escc <- escc[!grepl("N", escc$patientID), ]
escc_adj$patientID=gsub("N","",escc_adj$patientID)
escc$patientID=gsub("T","",escc$patientID)
wilcox.test(escc_adj$TLB_B_ratio, escc$TLB_B_ratio) #insignificant

clinic <- read.table(file="clinicInfo/ESCC_60patients.tsv",header=T)
table(clinic$Pathologic_stage)
table(clinic$N_stage)
table(clinic$M_stage)

escc <- merge(clinic, escc, by.x="SampleID", by.y="patientID")
escc$Pstage=as.factor(ifelse(escc$Pathologic_stage != "III", "Early", "Late"))
escc$Age_group=as.factor(ifelse(escc$Age<=60,"Young","Old"))
escc$Nstatus=as.factor(ifelse(escc$N_stage==0,"N","Y"))
for(i in c(2,4,5)){
    escc[,i]=as.factor(escc[,i])   
}
wilcox.test(escc$TLB_B_ratio ~ escc$Age_group)
#data:  escc$TLB_B_ratio by escc$Age_group
#W = 502, p-value = 0.03198
median(escc[escc$Age_group=="Young",]$TLB_B_ratio)
#0.008093525
median(escc[escc$Age_group=="Old",]$TLB_B_ratio)
#0.02285714

gej_se <- readRDS(file="data/GEJ_IM_TLB.rds")
gej=gej_se@meta.data
colnames(gej)[1]="patientID"
gej$patientID=gsub("-I","",gej$patientID)
gej$TLB=ifelse(gej$seurat_clusters==10,"Y","N")
count <- table(gej$patientID,gej$TLB)
gej <- data.frame(patientID=rownames(count),nTLB=unname(count[,2]), nOtherB=unname(count[,1]))
gej$TLB_B_ratio <- gej$nTLB/gej$nOtherB
gej_adj <- gej[grepl("N", gej$patientID), ]
gej <- gej[!grepl("N", gej$patientID), ]
gej_adj$patientID=gsub("N","",gej_adj$patientID)
gej$patientID=gsub("T","",gej$patientID)
wilcox.test(gej_adj$TLB_B_ratio, gej$TLB_B_ratio) #insignificant

clinic <- read.table(file="clinicInfo/GEJ_40patients.tsv",header=T)
table(clinic$Grade)
gej <- merge(clinic, gej, by.x="SampleID", by.y="patientID")
gej$Stage=as.factor(ifelse(gej$Grade<3,"Early", "Late"))
gej$Age_group=as.factor(ifelse(gej$Age<=60,"Young","Old"))
gej$Nstatus=as.factor(ifelse(!grepl("N0",gej$TNM),"Y","N"))
for(i in c(2,7,8)){
    gej[,i]=as.factor(gej[,i])
}
wilcox.test(gej$TLB_B_ratio ~ gej[,i])

df_combine <- rbind(escc[,c("SampleID","TLB_B_ratio","Age_group")],gej[,c("SampleID","TLB_B_ratio","Age_group")])
wilcox.test(df_combine$TLB_B_ratio ~ df_combine$Age_group)
#data:  df_combine$TLB_B_ratio by df_combine$Age_group
#W = 1237, p-value = 0.005417
median(df_combine[df_combine$Age_group=="Young",]$TLB_B_ratio)
# 0.007518797
median(df_combine[df_combine$Age_group=="Old",]$TLB_B_ratio)
# 0.01791486
