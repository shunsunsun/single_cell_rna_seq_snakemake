suppressPackageStartupMessages({
	library(infercnv)
})

args  <- commandArgs(trailingOnly=T)
type <- args[1] #escc, acgej, combined
workdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell"
indir <- paste0(workdir,"/CNV/",toupper(type),"_inferCNV")
outdir <- paste0(indir,"_v2")

objfile <- paste0(indir, "/infercnv_obj.rda")
if(file.exists(objfile)){
   load(file=objfile)
}else{
   infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(indir,"/count_matrix.txt"),
	annotations_file=paste0(indir, "/cell.annot.txt"),
	delim="\t", gene_order_file= paste0(indir, "/gene.annot.txt"),
	ref_group_names=c('ref-endo','ref-fib')) 
}

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=indir, 
                             cluster_by_groups=F, 
                             denoise=TRUE,
                             HMM=TRUE)

infercnv_obj2 = infercnv::run(infercnv_obj,
                             cutoff=0.1,  
                             out_dir=outdir) 
