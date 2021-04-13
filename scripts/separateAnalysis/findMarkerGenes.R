#Ref: https://hbctraining.github.io/scRNA-seq/lessons/09_merged_SC_marker_identification.html
#https://github.com/satijalab/seurat/issues/2180

suppressPackageStartupMessages({
	library(Seurat)
	library(future)
	library(tidyverse)
})

args  <- commandArgs(trailingOnly=T)
infile <- args[1]
sctransform <- args[2] #T or F
assay_used <- ifelse(sctransform=="T","SCT","RNA") #RNA is recommended; I only find SCT significantly different from RNA results and useful for integrated datasets
affix <- ifelse(sctransform=="T","_sct","")
outfile1 <- gsub(".rds",paste0(affix,"_conservedMarkers.rda"),infile)
outfile2 <- gsub(".rds",paste0(affix,"_allMarkers.rda"),infile)
#plotPrefix <- paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/cellident/",
#        gsub("_clustStab","", basename(dirname(infile))))
annotations <- read.csv("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/resources/annotation.csv")

nworker <- min(as.numeric(args[3]),length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3) #60G
randSeed <- 1129L
set.seed(randSeed)

se <- readRDS(file=infile)
clustRes <- colnames(se@meta.data)[grepl("_res\\.", colnames(se@meta.data))]
Idents(se) <- clustRes
#DefaultAssay(se) <- "RNA"

#Get conserved markers for any given cluster
#It is possible that when you run this function on all clusters, in some cases you will have clusters 
#that do not have enough cells for a particular group - and your function will fail. 
#For these clusters you will need to use FindAllMarkers()
get_conserved <- function(cluster){
    tryCatch({
	markers <- FindConservedMarkers(se, ident.1 = cluster,
        assay=assay_used, slot="data", grouping.var = "orig.ident", 
	only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	markers <- markers %>% rownames_to_column(var = "gene") %>%
	    left_join(y = unique(annotations[, c("gene_name", "description")]),
               by = c("gene" = "gene_name")) %>% cbind(cluster_id = cluster, .)
	return(markers)
    }, error=function(cond){
	print(paste0("Cluster", cluster, " has too few cells to run FindConservedMarkers"))
	return(NULL)
    })
}

conserved_markers_list <- list()
for(i in unique(se[[clustRes,drop=TRUE]])){
  markers_df <- get_conserved(i)
  if(!is.null(markers_df)){
	  conserved_markers_list[[paste0("Cluster", i)]] <- markers_df
	  print(x = head(markers_df))
	  #markers_genes = head(x = markers_df, n = 5)$gene
	  #VlnPlot(object = se,features=markers_genes,log=TRUE)
	  #ggsave(filename=paste0(plotPrefix,'_Cluster',i,'_markers_VlnPlot.pdf'))
	  #FeaturePlot(object = se, features=markers_genes, sort.cell=TRUE, min.cutoff='q10',label=TRUE, repel=TRUE)
	  #ggsave(filename=paste0(plotPrefix,'_Cluster',i,'_markers_FeaturePlot.pdf'))
  }
}
save(conserved_markers_list, file=outfile1)
all_markers <- FindAllMarkers(object = se, assay=assay_used, slot="data", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )
save(all_markers, file=outfile2)

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)


##Analyze the results
##Check analyzeMarkerGenes.R
