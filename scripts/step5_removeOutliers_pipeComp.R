input <- snakemake@input[[1]]
filename <- unname(unlist(strsplit(basename(input),"_")))
cancer <- filename[1]
celltype <- filename[2]
outdir=dirname(snakemake@output[["rds"]])
if(!dir.exists(outdir)){
        dir.create(outdir,recursive=T)
}

suppressPackageStartupMessages(library(Seurat))
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(scater))

## Read input
sce=readRDS(input)

####################
## Cell-level filtering

##customized filters used in isOutlier(): field_name:type:nmad
##percent_top_20 means the percent of counts in top-20 features
filters <- c("log10_total_counts:higher:2.5",
             "log10_total_counts:lower:5",
             "log10_total_features:higher:2.5",
             "log10_total_features:lower:5",
             "percent_top_20:both:5",
             "featcount_dist:both:5")

out <- lapply(strsplit(filters,":"), FUN=function(f){
	which(isOutlier(sce[[f[1]]], log=FALSE,
        	nmads=as.numeric(f[3]), type=f[2], batch=sce$orig.ident))
})

if(!snakemake@params[["all"]] && snakemake@params[["lowqual_list"]] != "NA"){
	lowqual_samples=read.table(file=snakemake@params[["lowqual_list"]],header=T)
	lowqual_samples=lowqual_samples[lowqual_samples$Cancer==cancer & lowqual_samples$Cell_type==celltype,]$Sample
	if(!is.null(lowqual_samples) && length(lowqual_samples)>0){
		subset=!(sce$orig.ident %in% lowqual_samples)
	}else{
		subset=NULL
	}
}else{		
	subset=NULL
}
mtout <- isOutlier(sce$mitoCountRatio, nmads=3, type="lower", subset=subset) |
	(isOutlier(sce$mitoCountRatio, nmads=2.5, type="higher", subset=subset) & 
	sce$mitoCountRatio > snakemake@params[["mtRatio"]])

out <- c(out, list(mt=which(mtout)))
out <- table(unlist(out))
out <- as.numeric(names(out)[which(out>=2)]) #cells deemed outliers in at least 2 of the thresholds
if(length(out)>0) sce <- sce[,-out]
sce=sce[Matrix::rowSums(counts(sce) > 0) >= snakemake@params[["nCell"]], Matrix::colSums(counts(sce) > 0) >= 0] 
	
filtered_seurat <- as.Seurat(sce)
saveRDS(filtered_seurat,file=snakemake@output[["rds"]])
