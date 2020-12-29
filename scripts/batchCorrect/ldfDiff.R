#Ref:https://www.bioconductor.org/packages/release/bioc/vignettes/CellMixS/inst/doc/CellMixS.html
#ldfDiff calculates the differences between each cell’s local density factor before and after data integration
#In an optimal case relative densities (according to the same set of cells) should not change by integration and the ldfDiff score should be close to 0. 
#In general the overall distribution of ldfDiff should be centered around 0 without long tails.
suppressPackageStartupMessages({
        library(CellMixS)
	library(SingleCellExperiment)
	library(cowplot)
	library(ggplot2)
	library(scater)
})
args  <- commandArgs(trailingOnly=T)
infile <- args[1]
#infile <- snakemake@input[[1]] #GEJ_QCed_*_walktrap_clustered.rds
plotDir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/batchCorrect/"
label <- args[2]
kVal <- as.numeric(args[3]) #70
useImmun=args[4]
set.seed(1129)

sce_noBC <- readRDS(file=paste0(dirname(infile),"/GEJ_QCed_sctNorm_noBatchCorr_walktrap_clustered.rds"))
sce <- readRDS(file=infile)
if(useImmun=='T'){
    meta <- colData(sce)
    sce <- sce[,grepl("-I",rownames(meta))]
    meta <- colData(sce_noBC)
    sce_noBC <- sce_noBC[,grepl("-I",rownames(meta))]
}
sce$orig.ident=as.factor(sapply(strsplit(sce$orig.ident,"-"),"[",1))
sce_noBC$orig.ident=sapply(strsplit(sce_noBC$orig.ident,"-"),"[",1)
split_sce_noBC <- list()
for(p in unique(sce_noBC$orig.ident)){
	split_sce_noBC[[p]] <- sce_noBC[,sce_noBC$orig.ident == p]
}
rm(sce_noBC)

for(r in reducedDimNames(sce)){
   if(length(reducedDimNames(sce))>1){
        label=paste0(label,"_",r)
        outfile <- gsub("walktrap_clustered",paste0(r,"_ldfK",kVal),infile)
   }else{
        outfile <- gsub("walktrap_clustered",paste0("ldfK",kVal),infile)
   }
   if(useImmun=='T'){
        outfile <- gsub("ldf","immCel_ldf",outfile)
   }
   if(!file.exists(outfile)){	
	ldf_res <- ldfDiff(sce_pre_list=split_sce_noBC, sce_combined=sce, k=kVal, group="orig.ident", res_name=label, dim_combined=r,
                dim_red="PCA", assay_pre="counts", n_dim=ncol(reducedDim(sce,r)))
	meta <- colData(ldf_res)
	ldf <- meta[,grepl("diff_ldf",colnames(meta))]
	res <- data.frame(BasedOn='all_cell',median_ldf=median(ldf),min_ldf=min(ldf),max_ldf=max(ldf))
	#saveRDS(ldf_res,file=gsub(".rds","_ldfResTmp.rds",outfile))
	
        walktrap_res <- which(grepl("snn",colnames(meta)) & !grepl("_res",colnames(meta)))
        colnames(meta)[walktrap_res]=gsub("snn","walktrap_snn",colnames(meta)[walktrap_res])
	colData(ldf_res) <- meta

	#plist <- list()
	for(i in which(grepl("snn", colnames(meta)) | grepl("_res.", colnames(meta)))){
                clust_res <- as.factor(meta[,i])
		clustMethod <- colnames(meta)[i]
		
		#plotfile2 <- paste0(plotDir, gsub(".rds",paste0("_",clustMethod,".jpg"),basename(outfile)))
		#print(plotfile2)
		#jpeg(plotfile2, type="cairo", width = 25, height=ceiling(length(unique(clust_res)))*20, units="cm", res=300)
		#plist[[clustMethod]] <- visIntegration(ldf_res, metric_var = paste0("diff_ldf.",label), cluster_var = clustMethod)
		#dev.off()		

                clust_ldf_median <- NULL
		clust_ldf_min <- NULL
		clust_ldf_max <- NULL
                for(c in unique(clust_res)){
                    clust_ldf <- meta[which(clust_res==c),grepl("diff_ldf",colnames(meta))]
                    clust_ldf_median <- c(clust_ldf_median, median(clust_ldf))
		    clust_ldf_min <- c(clust_ldf_min, min(clust_ldf))
		    clust_ldf_max <- c(clust_ldf_max, max(clust_ldf))
                }
                res <- rbind(res, data.frame(BasedOn=clustMethod,median_ldf=median(clust_ldf_median),
			min_ldf=median(clust_ldf_min),max_ldf=median(clust_ldf_max)))
	}
        saveRDS(res, file=outfile)
	#save(plist,file=paste0(plotDir, gsub(".rds","_plots.rda",basename(outfile))))
    }else{
        print(paste0("Already done: ", basename(outfile)))
    }
}
