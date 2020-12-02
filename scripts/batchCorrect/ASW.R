##The silhouette score of a data point is computed by subtracting its average distance to other members in the same cluster from its average distance to all members 
##of the neighboring clusters, and then dividing by the larger of the two values. The resulting score ranges from − 1 to 1, where a high score denotes that 
##the data point fits well in the current cluster, while a low score denotes a poor fit.

suppressPackageStartupMessages({
	library(Seurat)
	library(scater)
    	library(bluster)
	library(ggplot2)
	library(gridExtra)
})
args  <- commandArgs(trailingOnly=T)
infile=args[1]
#infile=snakemake@input[[1]] #GEJ_QCed_*_walktrap_clustered.rds
plotDir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/batchCorrect/"
sce <- readRDS(file=infile)
df <- as.data.frame(colData(sce))
df <- df[, colnames(df)=='orig.ident' | grepl("^snn", colnames(df)) | grepl("_res.", colnames(df))]
colnames(df) <- gsub("^snn","walktrap_snn",colnames(df))     
df$orig.ident=sapply(strsplit(df$orig.ident,"-"),"[",1)

for(r in reducedDimNames(sce)){
    outfile <- gsub("walktrap_clustered",paste0(r, "_asw"),infile)
    if(!file.exists(outfile)){
        rd <- reducedDim(sce, type=r)
	plist <- list()
	for(i in which(grepl("snn", colnames(df)) | grepl("_res.", colnames(df)))){
		clustMethod <- colnames(df)[i]
		sil_data <- as.data.frame(approxSilhouette(rd, clusters=df[,i]))
		sil_data$closest <- factor(ifelse(sil_data$width > 0, df[,i], sil_data$other))
           	sil_data$cluster <- factor(df[,i])
                res <- data.frame(BasedOn="all_cell", median=median(sil_data$width), min=min(sil_data$width), max=max(sil_data$width))
                clust_sil_median <- NULL
                clust_sil_min <- NULL
                clust_sil_max <- NULL
                for(c in unique(sil_data$cluster)){
                    clust_sil <- sil_data[sil_data$cluster==c, ]
                    clust_sil_median <- median(clust_sil$width)
                    clust_sil_min <- min(clust_sil$width)
                    clust_sil_max <- max(clust_sil$width)
                }
                res <- rbind(res, data.frame(BasedOn=clustMethod,median=median(clust_sil_median),min=median(clust_sil_min),max=median(clust_sil_max)))
		nClust <- length(unique(levels(sil_data$cluster)))
		if(nClust<=26){
			my_color_palette <- DiscretePalette(n=nClust,palette="alphabet")
		}else if(nClust<=36){ 
			my_color_palette <- DiscretePalette(n=nClust,palette="polychrome")
		}else{
			my_color_palette <- NULL
		}
		if(!is.null(my_color_palette)){
	                plist[[clustMethod]] <- ggplot(sil_data, aes(x=cluster, y=width, colour=closest)) +
        	                ggbeeswarm::geom_quasirandom(method="smiley",size=.5) + 
				scale_colour_manual(values = my_color_palette) + theme(legend.position="none") +
                	        labs(x=paste0("Clustering_",clustMethod),y="Approx. Silhouette width")
		}else{
			plist[[clustMethod]] <- ggplot(sil_data, aes(x=cluster, y=width, colour=closest)) +
				ggbeeswarm::geom_quasirandom(method="smiley",size=.5) + theme(legend.position="none") +
				labs(x=paste0("Clustering_",clustMethod),y="Approx. Silhouette width")
		}
	}
        saveRDS(res, file=outfile)
	plotfile <- paste0(plotDir, gsub("rds","jpg",basename(outfile)))
        if(!file.exists(plotfile)){
            if(length(plist)>3){
                jpeg(plotfile, width = 100, height = 60, units="cm", res=300)
                do.call("grid.arrange", c(plist, ncol=3))
                dev.off()
            }else{
                jpeg(plotfile, width = 60, height = 40, units="cm", res=300)
                do.call("grid.arrange", c(plist, ncol=2))
                dev.off()
            }
        }
    }else{
        print(paste0("Already done: ", basename(outfile)))
    }
}
