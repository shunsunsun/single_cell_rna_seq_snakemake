options(stringsAsFactors = F)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(reticulate))
if(!py_module_available(module = "scanorama")){
	warning("Cannot find python module scanorama. Make sure you have activated the right conda environment")
}else{
	scanorama <- import('scanorama')
	infile <- snakemake@input[[1]]
	outfile <- snakemake@output[[1]]

	# https://github.com/brianhie/scanorama/issues/38#issuecomment-551446738
	se <- readRDS(file=infile)
	se@meta.data$orig.ident=sapply(strsplit(se@meta.data$orig.ident,"-"),"[",1)
	assaylist <- list()
	genelist <- list()
	seuratobjectlist <- SplitObject(se,split.by = "orig.ident")
	for(i in 1:length(seuratobjectlist)){
		se <- seuratobjectlist[[i]]
		assaylist[[i]] <- t(as.matrix(GetAssayData(object=se[["SCT"]], slot="data")))
		genelist[[i]] <- rownames(se)
	}

	#integrated.data <- scanorama$integrate(assaylist, genelist)
	#corrected.data <- scanorama$correct(assaylist, genelist, return_dense=TRUE)
	integrated.corrected.data <- scanorama$correct(assaylist, genelist, return_dimred=TRUE, return_dense=TRUE)

	#In the returned list, element 1 is the dimensional reduction embeddings from Scanorama; 
	#element 2 is the corrected counts (element number 2); 
	#element 3 is common genes, which we use as rownames of the integrated, batch-corrected expression matrix.
	#First we need to transpose the matrices again in order to have them in the appropriate format

	intdata <- lapply(integrated.corrected.data[[2]], t)
	panorama <- do.call(cbind, intdata)
	rownames(panorama) <- as.character(integrated.corrected.data[[3]])
	colnames(panorama) <- unlist(sapply(assaylist, rownames))
	intdimred <- do.call(rbind, integrated.corrected.data[[1]])
	colnames(intdimred) <- paste0("PC_", 1:100)

	#We also add standard deviations in order to draw Elbow Plots in Seurat
	stdevs <- apply(intdimred, MARGIN = 2, FUN = sd)
	
	pan.seurat <- CreateSeuratObject(counts = panorama, assay = "pano")
  	#Adding metadata from all previous objects 
	pan.seurat@meta.data <- do.call(rbind, lapply(seuratobjectlist, function(x) x@meta.data))
    	# VERY IMPORTANT: make sure the rownames of your metadata slot are the same as the colnames of your integrated expression matrix 
	rownames(pan.seurat@meta.data) <- colnames(pan.seurat)
	rownames(intdimred) <- colnames(pan.seurat)
	pan.seurat[["pca"]] <- CreateDimReducObject(embeddings = intdimred, stdev = stdevs, key = "PC_", assay = "pano")

	saveRDS(pan.seurat, file=outfile)
}
