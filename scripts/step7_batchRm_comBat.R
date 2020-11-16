#need to run at source /appsnew/source/R.4.0.2share.sh
options(stringsAsFactors = F)
suppressPackageStartupMessages({
	library(sctransform)
	library(Seurat)
	library(sva)
})

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
#covariates <- snakemake@params[["covariates"]] #concatenated by "+"
covariates <- NULL
set.seed(1129)

se <- readRDS(file=infile)
se@meta.data$orig.ident=sapply(strsplit(se@meta.data$orig.ident,"-"),"[",1)
mod <- NULL
if(!is.null(covariates)){
    mod <- stats::model.matrix(stats::as.formula(paste0("~", covariates)), data = se@meta.data)
}
mat <- GetAssayData(object = se[["SCT"]], slot = "data") 

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

flted <- preprocess_big(mat)
row_max = matrixStats::rowMaxs(flted)
row_min = matrixStats::rowMins(flted)
if(sum(row_max==row_min)>0){
  flted = flted[row_max != row_min,]
}
rv_genes <- which((apply(flted, 1, var)!=0) == "FALSE")
rv_genes_names <- rownames(flted)[rv_genes]
flted <- flted[!(rownames(flted) %in% rv_genes_names),]

source("/gpfs2/gaog_pkuhpc/users/liny/transcriptome_analysis/tcga_gtex/script/runCombat.r")
my_iter=15
res_assay <- ComBatWrapper(data.frame(flted), se@meta.data$orig.ident, mod, my_iter)
saveRDS(res_assay, file=outfile)

#res_assay <- ComBat(dat = mat, batch = se@meta.data$orig.ident, mod = mod)
#saveRDS(CreateSeuratObject(raw.data = res_assay), file=outfile)
