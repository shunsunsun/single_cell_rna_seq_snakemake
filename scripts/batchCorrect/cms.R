#Ref:https://www.bioconductor.org/packages/release/bioc/vignettes/CellMixS/inst/doc/CellMixS.html
#cms scores are p.values from hypothesis testing, without any batch effect the p.value histogram should be flat. An increased number of very small p-values indicates the presence of a batch-specific bias within data.

suppressPackageStartupMessages({
        library(CellMixS)
	library(SingleCellExperiment)
	library(cowplot)
	library(ggplot2)
	library(scater)
        library(BiocParallel)
})
args  <- commandArgs(trailingOnly=T)
infile <- args[1]
#infile <- snakemake@input[[1]] #GEJ_QCed_*_walktrap_clustered.rds
plotDir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/batchCorrect/"
nworker <- as.numeric(args[2])
label <- args[3]
kVal <- as.numeric(args[4]) #30 or 50
useImmun=args[5]
set.seed(1129)

sce <- readRDS(file=infile)
if(useImmun=='T'){
    meta <- colData(sce)
    sce <- sce[,grepl("-I",rownames(meta))]
}
sce$orig.ident=sapply(strsplit(sce$orig.ident,"-"),"[",1)
nworker=min(nworker,parallel::detectCores())
print(paste0("Use ",nworker, " workers"))

for(r in reducedDimNames(sce)){

   if(length(reducedDimNames(sce))>1){
	label=paste0(label,"_",r)
	outfile <- gsub("walktrap_clustered",paste0(r,"_cmsK",kVal),infile)
   }else{
	outfile <- gsub("walktrap_clustered",paste0("cmsK",kVal),infile)
   }
   if(useImmun=='T'){
        outfile <- gsub("cms","immCel_cms",outfile)
   }
   if(!file.exists(outfile)){	
	
	cms_res <- cms(sce, k=kVal, group="orig.ident", res_name=label, dim_red=r, 
		n_dim=ncol(reducedDim(sce,r)), BPPARAM=MulticoreParam(nworker))
	meta <- colData(cms_res)
	cms_smooth <- meta[,grepl("cms_smooth",colnames(meta))]
	cms <- meta[,grepl("cms\\.",colnames(meta))]
	res <- data.frame(BasedOn='all_cell',median_cms=median(cms),median_smooth=median(cms_smooth),
                min_cms=min(cms),min_smooth=min(cms_smooth),max_cms=max(cms),max_smooth=max(cms_smooth))
	#saveRDS(cms_res,file=gsub(".rds","_cmsResTmp.rds",outfile))
	
	#plist <- list()
	#plotfile1 <- paste0(plotDir, gsub(".rds","_global.jpg",basename(outfile)))
	#jpeg(plotfile1, type="cairo", width = 25, height = 20, units="cm", res=300)
	#plist[["allHist"]] <- visHist(cms_res)
	#dev.off()

        walktrap_res <- which(grepl("snn",colnames(meta)) & !grepl("_res",colnames(meta)))
        colnames(meta)[walktrap_res]=gsub("snn","walktrap_snn",colnames(meta)[walktrap_res])
	colData(cms_res) <- meta

	for(i in which(grepl("snn", colnames(meta)) | grepl("_res.", colnames(meta)))){
                clust_res <- as.factor(meta[,i])
		clustMethod <- colnames(meta)[i]
		
		#plotfile2 <- paste0(plotDir, gsub(".rds",paste0("_",clustMethod,".jpg"),basename(outfile)))
		#print(plotfile2)
		#jpeg(plotfile2, type="cairo", width = 25, height=ceiling(length(unique(clust_res)))*20, units="cm", res=300)
		#plist[[clustMethod]] <- visCluster(cms_res, metric_var = paste0("cms.",label), cluster_var = clustMethod)
		#dev.off()		

                clust_cms_median <- NULL
		clust_cms_min <- NULL
		clust_cms_max <- NULL
		clust_cms_sm_median <- NULL
                clust_cms_sm_min <- NULL
                clust_cms_sm_max <- NULL
                for(c in unique(clust_res)){
                    clust_cms <- meta[which(clust_res==c),grepl("cms\\.",colnames(meta))]
		    clust_cms_sm <- meta[which(clust_res==c),grepl("cms_smooth\\.",colnames(meta))]
                    clust_cms_median <- c(clust_cms_median, median(clust_cms))
		    clust_cms_min <- c(clust_cms_min, min(clust_cms))
		    clust_cms_max <- c(clust_cms_max, max(clust_cms))
		    clust_cms_sm_median <- c(clust_cms_sm_median, median(clust_cms_sm))
                    clust_cms_sm_min <- c(clust_cms_sm_min, min(clust_cms_sm))
                    clust_cms_sm_max <- c(clust_cms_sm_max, max(clust_cms_sm))
                }
                res <- rbind(res, data.frame(BasedOn=clustMethod,median_cms=median(clust_cms_median),median_smooth=median(clust_cms_sm_median),
			min_cms=median(clust_cms_min), min_smooth=median(clust_cms_sm_min),max_cms=median(clust_cms_max),max_smooth=median(clust_cms_sm_max)))
	}
        saveRDS(res, file=outfile)
	#save(plist,file=paste0(plotDir, gsub(".rds","_plots.rda",basename(outfile))))
    }else{
        print(paste0("Already done: ", basename(outfile)))
    }
}
