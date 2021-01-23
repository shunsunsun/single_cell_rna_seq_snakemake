##The silhouette score of a data point is computed by subtracting its average distance to other members in the same cluster from its average distance to all members 
##of the neighboring clusters, and then dividing by the larger of the two values. The resulting score ranges from − 1 to 1, where a high score denotes that 
##the data point fits well in the current cluster, while a low score denotes a poor fit.

suppressPackageStartupMessages({
	library(Seurat)
	library(scater)
    	library(bluster)
	library(ggplot2)
	library(patchwork)
})
args  <- commandArgs(trailingOnly=T)
infile=args[1]
outfile=gsub(".rds","_asw.rds",infile)
plotfile=gsub("rds","jpg",outfile)
sce <- readRDS(file=infile)
sce <- as.SingleCellExperiment(sce)
meta <- colData(sce)

if(!file.exists(outfile)){
	rd <- reducedDim(sce, type="PCA")
	clustRes <- meta[,grepl("_res.", colnames(meta))]
	sil_data <- as.data.frame(approxSilhouette(rd, clusters=clustRes))
	sil_data$closest <- factor(ifelse(sil_data$width > 0, clustRes, sil_data$other))
        sil_data$cluster <- factor(clustRes)
        res <- data.frame(BasedOn="all_Cell", median=median(sil_data$width), min=min(sil_data$width), max=max(sil_data$width))
        for(c in unique(sil_data$cluster)){
            clust_sil <- sil_data[sil_data$cluster==c, ]
            clust_sil_median <- median(clust_sil$width)
            clust_sil_min <- min(clust_sil$width)
            clust_sil_max <- max(clust_sil$width)
            res <- rbind(res, data.frame(BasedOn=paste0("Cluster",c),median=clust_sil_median,min=clust_sil_min,max=clust_sil_max))
        }
        saveRDS(res, file=outfile)
	if(!file.exists(plotfile)){
		nClust <- length(unique(levels(sil_data$cluster)))
		if(nClust<=26){
			my_color_palette <- DiscretePalette(n=nClust,palette="alphabet")
		}else if(nClust<=36){ 
			my_color_palette <- DiscretePalette(n=nClust,palette="polychrome")
		}else{
			my_color_palette <- NULL
		}
		if(!is.null(my_color_palette)){
	      	        p <- ggplot(sil_data, aes(x=cluster, y=width, colour=closest)) +
        	        	        ggbeeswarm::geom_quasirandom(method="smiley",size=.5) + 
					scale_colour_manual(values = my_color_palette) + theme(legend.position="none") +
	                	        labs(x="Clusters",y="Approx. Silhouette width")
		}else{
			p <- ggplot(sil_data, aes(x=cluster, y=width, colour=closest)) +
					ggbeeswarm::geom_quasirandom(method="smiley",size=.5) + theme(legend.position="none") +
					labs(x="Clusters",y="Approx. Silhouette width")
			}
		}else{
			print("Skip plotting")
		}
		ggsave(filename=plotfile, plot=p, width=15, height=10, units="in", dpi=300)
	}
}else{
        print(paste0("Already done: ", basename(outfile)))
}
