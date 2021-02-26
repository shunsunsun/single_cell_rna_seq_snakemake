suppressPackageStartupMessages({
	library(Seurat)
	library(BiocParallel)
	library(ggplot2)
	library(patchwork)
})

samplefile <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/config/sample_cancertype.txt"
randseed=1129L

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
nworkers=min(as.numeric(args[2]),parallel::detectCores())
bpParam=BiocParallel::MulticoreParam(nworkers)

rerunUmap <- args[3]
nPC <- as.numeric(args[4])
k <- as.numeric(args[5])
reduc <- args[6]
object <- args[7] #se (seurat) or sce

if(grepl("Combined", infile)){
	type="combined"
}else if(grepl("ESCC",infile)){
	type="escc"
}else{
	type="gej"
}

prefix <- paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/cellident/",
       type, "_", gsub(".rds","_", basename(infile)))

##-------Get data------------
se <- readRDS(file=infile)
se$seurat_clusters <- NULL
if(object=="sce"){
	se <- as.Seurat(se)
}
if(!any(grepl("umap",names(se))) || rerunUmap == "T"){
	se <- RunUMAP(se, reduction=reduc, dims=1:nPC, n.neighbors=k)
	saveRDS(se,file=infile)
}	
meta <- se@meta.data
clustResCol <- colnames(meta)[grepl("_res\\.", colnames(meta))]
#my_color_palette <- DiscretePalette(n=length(unique(se[[clust, drop=TRUE]])), palette = "alphabet")
for(c in clustResCol){
	p1 = DimPlot(se, group.by=c, label=T, reduction="umap") + NoLegend()
	p2 = DimPlot(se, group.by="orig.ident", shuffle=T, seed=randseed, label=F, reduction="umap") + NoLegend()
	if(type=="combined"){
		gej_patients <- read.table(file=samplefile,header=T,stringsAsFactors=F)
		gej_patients <- unique(gej_patients[gej_patients$type=="G",]$patient)
		se[["cancer_type"]] <- ifelse(se[["orig.ident",drop=T]] %in% gej_patients, "G", "E")
		p3 = DimPlot(se, group.by="cancer_type", shuffle=T, seed=randseed, label=F, reduction="umap")
		p = (p1 + p2)/p3
		ggsave(filename=paste(prefix,c,".pdf"), plot=p, width = 10, height = 8, units="in",dpi=300)
	}else{
		p = p1 + p2
		ggsave(filename=paste(prefix,c,".pdf"), plot=p, width = 10, height = 4, units="in",dpi=300)
	}
}
