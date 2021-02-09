suppressPackageStartupMessages({
	library(biomaRt)
	library(dplyr)
	library(Seurat)
})

args  <- commandArgs(trailingOnly=T)
type <- args[1] #escc, acgej, combined
workdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell"

if(type=="escc"){
	indir <- paste0(workdir, "/data/ESCC_QCed_sctNorm_BatchHmy_clustStab")
	epi <- readRDS(paste0(indir,"/EpithelialCells.rds"))
	all <- readRDS(paste0(indir,"/Louvain_clust_k100_pc50_res0.2.rds"))
	endo_clust=c(7,16)
	fib_clust=c(2,6,13)
}else if(type=="acgej"){
	indir <- paste0(workdir,"/data/GEJ_QCed_sctNorm_BatchCCA_clustStab")
	epi <- readRDS(paste0(indir,"/EpithelialCells.rds"))
        all <- readRDS(paste0(indir,"/Louvain_clust_k100_pc50_res0.4.rds"))
        endo_clust=c(7)
        fib_clust=c(8,13)
}

outdir <- paste0(workdir,"/CNV/",toupper(type),"_inferCNV")
dir.create(outdir, showWarnings=F)

epi.source <- ifelse(grepl("N", epi@meta.data$orig.ident), "adj_epi", "tumor_epi")
epiMat <- as.data.frame(GetAssayData(epi,assay="RNA",slot='counts'))
rm(epi)

meta <- all@meta.data
clustCol <- which(grepl("_res\\.",colnames(meta)))
colnames(meta)[clustCol]="clustRes"
fib <- rownames(meta[meta$clustRes %in% fib_clust, ])
endo <- rownames(meta[meta$clustRes %in% endo_clust, ])

fib=sample(fib,800)
endo=sample(endo,800)
fibMat=as.data.frame(GetAssayData(subset(all,cells=fib), assay="RNA", slot='counts'))
endoMat=as.data.frame(GetAssayData(subset(all,cells=endo), assay="RNA", slot='counts'))

dat=cbind(epiMat,fibMat,endoMat)
cell.annot=data.frame(v1=colnames(dat),v2=c(epi.source,
	rep('spike-fib',300),rep('ref-fib',500),rep('spike-endo',300), rep('ref-endo',500)))
write.table(cell.annot, file=paste0(outdir, "/cell.annot.txt"), row.names=F, col.names = F, quote=F, sep = "\t" )

query_genes = rownames(dat)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
results_id <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", 'chromosome_name',
      'start_position', 'end_position'),filters = "hgnc_symbol", values=query_genes, mart = ensembl)

chromo_list <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                 "11", "12", "13", "14", "15", "16", "17", "18", "19",
                 "20", "21", "22", "X", "Y")

results <- results_id %>% filter(chromosome_name %in% chromo_list)
results <- results %>% select(hgnc_symbol, chromosome_name, start_position, end_position)

#check if any duplicates in gene position
rep_gene <- data.frame(table(results$hgnc_symbol))
results[results$hgnc_symbol %in% rep_gene[rep_gene$Freq>1, ]$Var1, ]

#clear replicates
results_unique <- results[!duplicated(results$hgnc_symbol), ]

# write table of gene notations
write.table(results_unique, file=paste0(outdir, "/gene.annot.txt"), col.names = FALSE, row.names = FALSE,  sep = "\t", quote = FALSE)

# filter the counts matrix according to results of chromosome positions
dat <- dat[rownames(dat) %in% results_unique$hgnc_symbol,]
dat <- dat[match(results_unique$hgnc_symbol, rownames(dat)),] 
dim(dat)

write.table(dat,file=paste0(outdir,"/count_matrix.txt"),sep="\t",quote=F)
#saveRDS(dat, file = paste0(outdir,"/count_matrix.rds"))
