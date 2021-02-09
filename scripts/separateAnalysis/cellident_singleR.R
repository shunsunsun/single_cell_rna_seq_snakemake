suppressPackageStartupMessages({
	library(Seurat)
	library(SingleR)
	library(celldex)
	library(BiocParallel)
	library(ggplot2)
	library(patchwork)
})

ref <- "hpca_be_dice_nh_mi"
label <- "label.main" #label.main or label.fine
refdir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/resources/SingleRdata"

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
nworkers=min(as.numeric(args[2]),ifelse(grepl("ESCC",basename(infile)), 10, parallel::detectCores()))
bpParam=BiocParallel::MulticoreParam(nworkers)

rerunUmap <- args[3]
nPC <- as.numeric(args[4])
k <- as.numeric(args[5])
reduc <- args[6]

outfile <- gsub(".rds",paste0("_singleR_",label,".rda"),infile)
plotfile <- paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/cellident/",
       gsub("clustStab","", basename(dirname(infile))), "PC", nPC, "_k", k, ".pdf")
#	gsub("clustStab","", basename(dirname(infile))), "PC", nPC, "_k", k, "_singleR.pdf") 
#tmpfile <- gsub("pdf","rda",plotfile)

#if(!file.exists(plotfile)){
	if(!file.exists(outfile)){
		##--------Get reference---------

		ref_list <- list()
		label_list <- list()
		if(grepl("hpca",ref)){
			ref_list[["hpca"]] <- readRDS(file=paste0(refdir,'/', 'HumanPrimaryCellAtlasData.rds'))
			label_list[["hpca"]] <- ref_list[["hpca"]][[label]]
		}
		if(grepl("be",ref)){
			ref_list[["be"]] <- readRDS(file=paste0(refdir,'/', 'BlueprintEncodeData.rds'))
			label_list[["be"]] <- ref_list[["be"]][[label]]
		}
		if(grepl("dice",ref)){
			ref_list[["dice"]] <- readRDS(file=paste0(refdir,'/', 'DatabaseImmuneCellExpressionData.rds'))
			label_list[["dice"]] <- ref_list[["dice"]][[label]]
		}
		if(grepl("nh",ref)){
			ref_list[["nh"]] <- readRDS(file=paste0(refdir,'/', 'NovershternHematopoieticData.rds'))
			label_list[["nh"]] <- ref_list[["nh"]][[label]]
		}
		if(grepl("mi",ref)){
			ref_list[["mi"]] <- readRDS(file=paste0(refdir,'/', 'MonacoImmuneData.rds'))
			label_list[["mi"]] <- ref_list[["mi"]][[label]]
		}

		##-------Get data------------
		se <- readRDS(file=infile)
		meta <- se@meta.data
		test=GetAssayData(object = se[["SCT"]], slot = "data")
		clustRes <- meta[,grepl("_res\\.",colnames(meta))]
		rm(meta)
		if(length(clustRes)>0){
			##------SingleR identification--------	
			#Annotation are performed on the cluster-level profiles rather than on the single-cell level. it directly returns the likely cell type identity of each cluster
			#This approach assumes that each cluster in the test dataset corresponds to exactly one reference label, which is likely not true
			clust.pred <- SingleR(test = test, ref = ref_list, labels = label_list, recompute=TRUE, clusters=clustRes, BPPARAM=bpParam)
			save(clust.pred, file=outfile)
		}else{
			print("Error: no clustering info in the metadata")
		}
	}else{
		load(file=outfile)
		se=readRDS(file=infile)
	}
	if(!any(grepl("umap",names(se))) || rerunUmap == "T"){
		se <- RunUMAP(se, reduction=reduc, dims=1:nPC, n.neighbors=k)
		saveRDS(se,file=infile)
	}	
	meta <- se@meta.data
	clustRes <- colnames(meta)[grepl("_res\\.", colnames(meta))]
        se[["SingleR_label"]] <- clust.pred$labels[match(clustRes, rownames(clust.pred))]
	se[["facs"]] <- ifelse(grepl("-E", meta$ident), "CD45-", "CD45+")
	#my_color_palette <- DiscretePalette(n=length(unique(se[[clust, drop=TRUE]])), palette = "alphabet")
	p1 = DimPlot(se, group.by=clustRes, label=T, reduction="umap") + NoLegend()
	p2 = DimPlot(se, group.by="SingleR_label", label=T, label.size=3, reduction="umap") + NoLegend()
#	ggsave(filename=plotfile, plot=p2, width=5, height=4, units="in", dpi=300)
	p3 = DimPlot(se, group.by="facs", label=F, reduction="umap")
	p4 = FeaturePlot(se,reduction="umap",features='PTPRC',order=T,min.cutoff='q10',repel=T)
	#save(p1, p2, p3, p4, file=tmpfile)
	p = (p1 + p3) / (p2 + p4)
   	p = p + plot_annotation(tag_levels = 'A') 
		#& theme(legend.text=element_text(size=8))
   	ggsave(filename=plotfile, plot=p, width = 10, height = 8, units="in",dpi=300)
}
