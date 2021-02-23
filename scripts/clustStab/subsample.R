suppressPackageStartupMessages({
	library(sctransform)
        library(Seurat)
	library(tidyverse)
	library(future)
})

## see https://bitbucket.org/snakemake/snakemake/issues/917/enable-stdout-and-stderr-redirection
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

subset_seurat_obj<- readRDS(snakemake@input[[1]])
k<- snakemake@wildcards[["k"]]
pc.use<- snakemake@wildcards[["pc"]]
run_id<- snakemake@wildcards[["run_id"]]
rate <- snakemake@params[["rate"]]
PreprocessSubsetData_pars<- snakemake@params[["PreprocessSubsetData_pars"]]
resString <- as.character(snakemake@params[["res"]])
resVal <- as.numeric(unname(unlist(strsplit(resString,','))))
if(length(resVal)>1){
        resString=paste0("c(",resString,")")
}
reducString <- snakemake@params[["reduc"]]
nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
output <- snakemake@output[[1]]
rseed=1129L

if(snakemake@params[["regressNum"]]){
    regVars=c('mitoCountRatio', 'nFeature_RNA', 'nCount_RNA')
}else{
    regVars=c('mitoCountRatio')
}

RandomSubsetData<- function(object, rate, random.subset.seed = NULL, ...){
        ncells<- nrow(object@meta.data)
        ncells.subsample<- round(ncells * rate)

        set.seed(random.subset.seed)

        selected.cells<- sample(colnames(object), ncells.subsample)
        object<- subset(object, cells =  selected.cells,...)
	DefaultAssay(object) <- "RNA"
	object <- DietSeurat(object, assays="RNA")
        return(object)
}

PreprocessSubsetData<- function(object,
                                variable.features.n = 3000,
                                num.pc = 50,
                                pc.use = NULL,
                                score.thresh = 1e-5,
                                sig.pc.thresh = 0.05,
                                n.start = 100,
                                nn.eps = 0,
                                resolution = seq(0.4,0.8,by=0.2),
				#resolution = 0.4,
				reduc = "pca",
                                k.param = 30,
				nworker = 8,
				random.seed = 1129L,
                                ...){

	
	cat(sprintf("Use %d workders with random seed %d\n",nworker,random.seed))
	plan("multiprocess", workers = nworker)
	options(future.globals.maxSize = 60*1024^3)

	if(!grepl("glmpca", reduc)){
	        regVars <- intersect(regVars, colnames(object[[]]))
        	# SCTransform replaces NormalizeData, ScaleData and FindVariableFeatures
	        object<- SCTransform(object, vars.to.regress = regVars, verbose = FALSE)

		if(reduc == "pca"){
		        object<- RunPCA(object = object, features = VariableFeatures(object = object), npcs = num.pc)
		}
		if(reduc == "harmony"){
			suppressPackageStartupMessages(library(harmony))
			object@meta.data$orig.ident=sapply(strsplit(object@meta.data$orig.ident,"-"),"[",1)
			object <- RunPCA(object = object, verbose=F)		
			object <- RunHarmony(object=object, assay.use = "SCT", dims.use = 1:50, reduction = "pca", group.by.vars = "orig.ident", plot_convergence = FALSE)	
		}
	}else{
		suppressPackageStartupMessages({
		        library(SeuratWrappers)
		        library(glmpca)
		})
		reduc <- paste0(reduc,pc.use)			
		m <- GetAssayData(object, slot = "counts", assay = "RNA")
		devs <- scry::devianceFeatureSelection(m)
		dev_ranked_genes <- rownames(object)[order(devs, decreasing = TRUE)]
		ndim <- pc.use
	        #Sparse matrices are coerced to dense matrice for  minibatch='none'; If this exhausts memory, consider setting minibatch to 'stochastic' or 'memoized'
        	object <- RunGLMPCA(object, features = head(dev_ranked_genes,n=2000), L = ndim, minibatch='stochastic',
	                reduction.name=reduc, reduction.key=paste0(toupper(reduc),"_"))
	}

        object <- FindNeighbors(object, dims = 1:pc.use, k.param = k.param, nn.eps = nn.eps, verbose = FALSE, reduction = reduc, force.recalc = TRUE)
        object <- FindClusters(object = object, n.start = n.start, resolution = resolution, random.seed=random.seed, method="igraph", verbose = FALSE)
	
	options(future.globals.maxSize = 500*1024^2) #500M
	plan(sequential)
	
        return(object)
}


subset_seurat_obj <- RandomSubsetData(subset_seurat_obj, rate = rate)

full_sample_clust <- list()
for (c in colnames(subset_seurat_obj@meta.data)[grepl("_res.",colnames(subset_seurat_obj@meta.data))]){
	r <- unname(unlist(strsplit(c,'res\\.')))[2]
	full_sample_clust[[ paste0("res",r) ]] <- subset_seurat_obj[[c]]
	subset_seurat_obj[[c]] <- NULL #remove old clustering results to get ready for new ones
}
if("seurat_clusters" %in% colnames(subset_seurat_obj@meta.data)){
	subset_seurat_obj$seurat_clusters <- NULL
}

command<- paste("PreprocessSubsetData", "(", "subset_seurat_obj,", "k.param=", k, ",", "pc.use=", pc.use, ",", "resolution=", resString, ",",
	"nworker=", nworker, ",", "random.seed=", rseed, ",", reducString, ",", PreprocessSubsetData_pars, ")")
print(paste0("Command=",command))
subset_seurat_obj<- eval(parse(text=command))

tmp_list <- list()
i=1
for (c in colnames(subset_seurat_obj[[]])[grepl("_res.",colnames(subset_seurat_obj[[]]))]){
    r <- unname(unlist(strsplit(c,'res\\.')))[2]
    tmp_list[[i]] <- tibble::tibble(pc = pc.use, resolution = r, k_param = k, original_ident_full = list(full_sample_clust[[ paste0("res",r) ]]),
	recluster_ident = list(subset_seurat_obj[[c]]), round = run_id)
    i <- i+1
}
res <- do.call(bind_rows, tmp_list)
saveRDS(res, file = output)

## make sure it is not empty file
#info<- file.info(output)
#if (info$size == 0) {
#    quit(status = 1)
#}
