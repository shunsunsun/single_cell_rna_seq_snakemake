##https://www.cellphonedb.org/explore-sc-rna-seq

library(Seurat)
library("biomaRt")
library(data.table)

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
outprefix <- args[2]
assay <- args[3] #SCT or RNA

se <- readRDS(infile)
counts <- GetAssayData(se, assay=assay, slot="data") #normalized counts
meta <- se@meta.data
rm(se)

###### break sparse matrix into trunk, convert into dense matrix and then concatenate
preprocess_big <- function(sparse_matrix,by=10000){
  n_col = ncol(sparse_matrix)
  n = floor(n_col/by)
  if(n<1){
    return(as.matrix(sparse_matrix))
  }else{
    for(i in 1:n){
      mat <- as.matrix(sparse_matrix[,((i-1)*by+1):(i*by)])      
      if(i<2){
        res <- mat
      }else{
        res <- cbind(res,mat)
      }
      rm(mat)
    }
    if(n_col > n*by){
      res <- cbind(res,as.matrix(sparse_matrix[,(n*by+1):n_col]))
    }
  }
  rm(sparse_matrix)
  return(res)
}

counts <- preprocess_big(counts)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
genes  <-  getBM(filters='hgnc_symbol',attributes = c('ensembl_gene_id','hgnc_symbol'),
	values = rownames(counts),mart = ensembl)
counts <- counts[rownames(counts) %in% genes$hgnc_symbol,]
counts <- tibble::rownames_to_column(as.data.frame(counts), var = 'hgnc_symbol')
counts <- merge(genes, counts, by="hgnc_symbol")
counts <- counts[,-1]
colnames(counts)[1]="Gene"
fwrite(counts, paste0(outprefix,'_count.txt'), sep='\t', quote=F, row.names=F, col.names=T)

clustRes <- paste0("C", meta[,grepl("_res\\.",colnames(meta))])
meta_data <- cbind(cell=rownames(meta), cell_type=clustRes)   #####  cluster is the user’s specific cluster column
fwrite(meta_data, paste0(outprefix,'_meta.txt'), sep='\t', quote=F, row.names=F, col.names=T)

