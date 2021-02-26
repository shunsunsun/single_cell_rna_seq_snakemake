#Reference list (in /gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/resources/scibet_scRNAseq_models)
#1 bloodcells_GSE94820.csv
#2 breastcancer_epithelialCells_GSE113197.csv
#3 breastcancer_tme_GSE75688.csv
#4 CD127positive_innate_lymphoid_cells.csv
#5 dendriticCells_GSE89232.csv
#6 dermal_fibroblasts_GSE104225.csv
#7 epithelialCells_GSE86618.csv
#8 lympho_myeloid_progeniter_GSE100618.csv
#9 major_human_cell_types.csv
#10 melanoma_cancercell_GSE115978.csv
#11 melanoma_non-cancercell_GSE115978.csv
#12 plasmaCells_GSE110499.csv
#13 prostatecancer_epithelial_GSE99795.csv
#14 stromal_B_NK_GSE72056.csv


suppressPackageStartupMessages({
	library(Seurat)
	library(ggplot2)
	library(tidyverse)
	library(scibet)
	library(viridis)
	library(ggsci)
})
#options(scibet.n_threads=10)
args  <- commandArgs(trailingOnly=T)
infile <- args[1]
refidx <- args[2] #use the list above
outfile <- args[3]
#outfile <- gsub(".rds","_scibet.rda",infile)
refdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/resources/scibet_scRNAseq_models/"
refs <- list.files(path=refdir, pattern="*.csv")
refs <- refs[as.numeric(unlist(strsplit(refidx,',')))]


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

se <- readRDS(file=infile)
query.expr <- GetAssayData(object = se, slot = "data", assay="SCT")
query.expr <- preprocess_big(query.expr)
rm(se)
query.expr <- t(query.expr)

ci <- list()
for(f in refs){
	ref_label=gsub(".csv","",f)
	model <- readr::read_csv(paste0(refdir,f))
	if(colnames(model)[1] == "X1" & colnames(model)[2] == "Unnamed: 0"){
     		model <- model[,-1]
    	} 
	model <- pro.core(model)
        prd <- LoadModel(model)
	ci[[ref_label]] <- prd(query.expr)
}
save(ci,file=outfile)
