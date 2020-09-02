library(Seurat)

input <- snakemake@input[[1]]
label <- gsub("_merged.rds","",basename(input))
outdir <- snakemake@output[[1]]

## Read input
seurat_obj <- readRDS(input)

seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
seurat_obj$mitoCountRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")/100
seurat_obj$riboCountRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^RP[SL][[:digit:]]")/100

library(scater)
sce=as.SingleCellExperiment(seurat_obj)
is.mito <- grepl("^MT-",rownames(sce))
sce=addPerCellQC(sce,subsets=list(Mito=is.mito))
#df=perCellQCMetrics(sce,subsets=list(Mito=is.mito))
#colData(sce)=cbind(colData(sce),df)

if(snakemake@params[["all_sample"]]=="true"){
        qc.lib <- isOutlier(sce$sum, log=TRUE, type="lower", batch=sce$orig.ident)
        qc.nexpr <- isOutlier(sce$detected, log=TRUE, type="lower", batch=sce$orig.ident)
        qc.mito <- isOutlier(sce$subsets_Mito_percent, type="higher", batch=sce$orig.ident)
} else if(snakemake@params[["lowqual_sample_list"]]!="NA"){
        qc.lib <- isOutlier(sce$sum, log=TRUE, type="lower", batch=sce$orig.ident)
        qc.nexpr <- isOutlier(sce$detected, log=TRUE, type="lower", batch=sce$orig.ident)
        qc.mito <- isOutlier(sce$subsets_Mito_percent, type="higher", batch=sce$orig.ident,
                subset=!(sce$orig.ident %in% lowqual_samples))
}
sce$discard_mad = qc.lib | qc.nexpr | qc.mito

if(snakemake@params[['outlyingness']]=="true"){
        library(robustbase)
        stats <- cbind(log10(sce$sum), log10(scef$detected), sce$subsets_Mito_percent)
        outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
        if(snakemake@params[["all_sample"]]=="true"){
                sce$discard_outlying <- isOutlier(outlying, type = "higher", batch=sce$orig.ident)
        }else if(snakemake@params[["lowqual_sample_list"]]!="NA"){
                sce$discard_outlying <- isOutlier(outlying, type = "higher", batch=sce$orig.ident,
                        subset=!(sce$orig.ident %in% lowqual_samples))
        }
}

sce_col_names=colnames(colData(sce))
for(i in which(grepl("discard",sce_col_names))){
	jpeg(sprintf("%s/%s_scater_outlier_nCount.jpg", outdir, label), width = 8, height = 6, units="in", res=300) 
	p=plotColData(sce, x="orig.ident", y="sum", colour_by=sce_col_names[i]) + xlab("") + ylab("log_sum") + 
		theme(axis.text.x=element_text(angle=90)) + scale_y_log10() + ggtitle("Total count per cell")
	print(p)
	dev.off()

	jpeg(sprintf("%s/%s_scater_outlier_nFeature.jpg", outdir, label), width = 8, height = 6, units="in", res=300) 
	p=plotColData(sce, x="orig.ident", y="detected", colour_by=sce_col_names[i]) + xlab("") + ylab("log_detected") + 
		theme(axis.text.x=element_text(angle=90)) + scale_y_log10() + ggtitle("Number of expressed genes per cell")
	print(p)
	dev.off()

	jpeg(sprintf("%s/%s_scater_outlier_mitoPerc.jpg", outdir, label), width = 8, height = 6, units="in", res=300)
	p=plotColData(sce, x="orig.ident", y="subsets_Mito_percent", colour_by=sce_col_names[i]) + xlab("") + 
		theme(axis.text.x=element_text(angle=90)) + ggtitle("Mito percent")
	print(p)
	dev.off()

	jpeg(sprintf("%s/%s_scater_outlier_nCount_mitoPerc.jpg", outdir, label), width = 8, height = 6, units="in", res=300)
	p=plotColData(sce, x="sum", y="subsets_Mito_percent", colour_by=sce_col_names[i]) + xlab("log_sum") +
		scale_x_log10() + ggtitle("Mito percent vs. total count")
	print(p)
	dev.off()

	jpeg(sprintf("%s/%s_scater_outlier_nCount_mitoPerc_perSample.jpg", outdir, label), width = 15, height = 5, units="in", res=300)
	p=ggplot(data=colData(sce), aes(x=sum, y=subsets_Mito_percent, color=sce_col_names[i]) + geom_point(size=0.1) + scale_x_log10() +
		stat_smooth() + theme_classic() + facet_wrap(~orig.ident,ncol=9,scales="free")
	print(p)
	dev.off()
}
	

	
