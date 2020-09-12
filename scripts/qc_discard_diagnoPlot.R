input <- snakemake@input[[1]]
filename <- unname(unlist(strsplit(basename(input),"_")))
cancer <- filename[1]
celltype <- filename[2]
label <- paste0(cancer,"_",celltype)
outdir <- snakemake@params[["outdir"]]
if(!dir.exists(outdir)){
	dir.create(outdir,recursive=T)
}


library(Seurat)
options(stringsAsFactors=F)

## Read input
seurat_obj <- readRDS(input)

library(scater)
sce=as.SingleCellExperiment(seurat_obj)
is.mito <- grepl("^MT-",rownames(sce))
sce=addPerCellQC(sce,subsets=list(Mito=is.mito))
#df=perCellQCMetrics(sce,subsets=list(Mito=is.mito))
#colData(sce)=cbind(colData(sce),df)

if(snakemake@params[["all"]]){
	strategy="_discard_allbased"
        qc.lib <- isOutlier(sce$sum, log=TRUE, type="lower", batch=sce$orig.ident)
        qc.nexpr <- isOutlier(sce$detected, log=TRUE, type="lower", batch=sce$orig.ident)
        qc.mito <- isOutlier(sce$subsets_Mito_percent, type="higher", batch=sce$orig.ident)
	sce$discard = qc.lib | qc.nexpr | qc.mito
} else if(snakemake@params[["lowqual_list"]] != "NA"){
	strategy="_discard_hqbased"
        qc.lib <- isOutlier(sce$sum, log=TRUE, type="lower", batch=sce$orig.ident)
        qc.nexpr <- isOutlier(sce$detected, log=TRUE, type="lower", batch=sce$orig.ident)
	lowqual_samples=read.table(file=snakemake@params[["lowqual_list"]],header=T)
	lowqual_samples=lowqual_samples[lowqual_samples$Cancer==cancer & lowqual_samples$Cell_type==celltype,]$Sample
	if(!is.null(lowqual_samples) && length(lowqual_samples)>0){
		qc.mito <- isOutlier(sce$subsets_Mito_percent, type="higher", batch=sce$orig.ident,
		subset=!(sce$orig.ident %in% lowqual_samples))
	}else{
        	qc.mito <- isOutlier(sce$subsets_Mito_percent, type="higher", batch=sce$orig.ident)
	}
	sce$discard = qc.lib | qc.nexpr | qc.mito
}

##We cannot do the following by excluding some samples
if(snakemake@params[["all"]] && snakemake@params[["outlyingness"]]){
        library(robustbase)
        stats <- cbind(log10(sce$sum), log10(scef$detected), sce$subsets_Mito_percent)
        outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
        sce$discard_outlying <- isOutlier(outlying, type = "higher", batch=sce$orig.ident)
}

sce_col_names=colnames(colData(sce))
for(i in which(grepl("discard",sce_col_names))){
	jpeg(sprintf("%s/%s%s_nCount.jpg", outdir, label, strategy), width = 8, height = 6, units="in", res=300) 
	p=plotColData(sce, x="orig.ident", y="sum", colour_by=sce_col_names[i]) + xlab("") + ylab("nUMI") +
		theme(axis.text.x=element_text(angle=90)) + scale_y_log10()
	print(p)
	dev.off()

	jpeg(sprintf("%s/%s%s_nFeature.jpg", outdir, label, strategy), width = 8, height = 6, units="in", res=300) 
	p=plotColData(sce, x="orig.ident", y="detected", colour_by=sce_col_names[i]) + xlab("") + ylab("nGene") + 
		theme(axis.text.x=element_text(angle=90)) + scale_y_log10()
	print(p)
	dev.off()

	jpeg(sprintf("%s/%s%s_mitoPerc.jpg", outdir, label, strategy), width = 8, height = 6, units="in", res=300)
	p=plotColData(sce, x="orig.ident", y="subsets_Mito_percent", colour_by=sce_col_names[i]) + xlab("") + ylab("mito percent") + 
		theme(axis.text.x=element_text(angle=90))
	print(p)
	dev.off()

	jpeg(sprintf("%s/%s%s_nCount_mitoPerc.jpg", outdir, label, strategy), width = 8, height = 6, units="in", res=300)
	p=plotColData(sce, x="sum", y="subsets_Mito_percent", colour_by=sce_col_names[i]) + xlab("nUMI") + ylab("mito percent") + 
		scale_x_log10()
	print(p)
	dev.off()

	jpeg(sprintf("%s/%s%s_nCount_mitoPerc_perSample.jpg", outdir, label, strategy), width = 30, height = 10, units="in", res=300)
	select_col=sym(sce_col_names[i])
	p=ggplot(data=as.data.frame(colData(sce)), aes(x=sum, y=subsets_Mito_percent, color=!!select_col)) + geom_point(size=0.1) + scale_x_log10() +
		xlab("nUMI") + ylab("mito percent") + theme(legend.position="bottom") + facet_wrap(~orig.ident,ncol=9,scales="free")
	print(p)
	dev.off()
}
	

	
