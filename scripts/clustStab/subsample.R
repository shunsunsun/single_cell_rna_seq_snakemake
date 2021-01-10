suppressPackageStartupMessages({
	library(Seurat)
	library(tidyverse)
	library(future)
})

## see https://bitbucket.org/snakemake/snakemake/issues/917/enable-stdout-and-stderr-redirection
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

seurat_obj<- readRDS(snakemake@input[[1]])
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
nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
output <- snakemake@output[[1]]
rseed=1129L

#if(snakemake@params[["regressNum"]]){
#    regVars=c('mitoCountRatio', 'nFeature_RNA', 'nCount_RNA')
#}else{
#    regVars=c('mitoCountRatio')
#}


RandomSubsetData<- function(object, rate, random.subset.seed = NULL, ...){
        ncells<- nrow(object@meta.data)
        ncells.subsample<- round(ncells * rate)

        set.seed(random.subset.seed)

        selected.cells<- sample(colnames(object), ncells.subsample)
        object<- subset(object, cells =  selected.cells,...)
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
                                k.param = 30,
				nworker = 8,
				random.seed = 1129L,
                                ...){

	
	cat(sprintf("Use %d workders with random seed %d\n",nworker,random.seed))
	plan("multiprocess", workers = nworker)
	options(future.globals.maxSize = 60*1024^3)

        #meta.data.colnames<- object@meta.data %>% colnames()
        #vars.to.regress<- c("percent.mt","nFeature_RNA")
        # in case the seurat object does not have percent.mito in metadata
        #vars.to.regress<- vars.to.regress[vars.to.regress %in% meta.data.colnames]
        # default is on variable features only, omit the features argument
        # SCTransform replaces NormalizeData, ScaleData and FindVariableFeatures
        #object<- SCTransform(object, vars.to.regress = vars.to.regress, variable.features.n = variable.features.n, verbose = FALSE)

        object<- RunPCA(object = object, features = VariableFeatures(object = object), npcs = num.pc)

        if (is.null(pc.use)){
                object<- JackStraw( object = object, num.replicate = 100, dims = num.pc)

                object <- ScoreJackStraw(object = object, dims = 1:num.pc, score.thresh = score.thresh)

                PC_pvalues<- object@reductions$pca@jackstraw@overall.p.values

                ## determin how many PCs to use.
                pc.use<- min(which(PC_pvalues[,"Score"] > sig.pc.thresh)) -1

        }

        # add significant pc number to metadata, need to have names same as the cells
        #pc.use.meta<- rep(pc.use, length(colnames(object)))
        #names(pc.use.meta)<- colnames(object)
        #object<- AddMetaData(object = object, metadata = pc.use.meta, col.name = "pc.use")
        object<- FindNeighbors(object, dims = 1:pc.use, k.param = k.param, nn.eps = nn.eps, verbose = FALSE, reduction = "pca", force.recalc = TRUE)
        object <- FindClusters(object = object, reduction = "pca", n.start = n.start, resolution = resolution, 
		random.seed=random.seed, method="igraph", verbose = FALSE)
	
	options(future.globals.maxSize = 500*1024^2) #500M
	plan(sequential)
	
        return(object)
}


subset_seurat_obj <- RandomSubsetData(seurat_obj, rate = rate)

full_sample_clust <- list()
for (c in colnames(subset_seurat_obj@meta.data)[grepl("_res.",colnames(subset_seurat_obj@meta.data))]){
	full_sample_clust[[c]] <- subset_seurat_obj[[c]]
	subset_seurat_obj[[c]] <- NULL #remove old clustering results to get ready for new ones
}
subset_seurat_obj$seurat_clusters <- NULL

command<- paste("PreprocessSubsetData", "(", "subset_seurat_obj,", "k.param=", k, ",", "pc.use=", pc.use, ",", "resolution=", resString, ",",
	"nworker=", nworker, ",", "random.seed=", rseed, ",", PreprocessSubsetData_pars, ")")
subset_seurat_obj<- eval(parse(text=command))

tmp_list <- list()
i=1
for (c in colnames(subset_seurat_obj[[]])[grepl("_res.",colnames(subset_seurat_obj[[]]))]){
    r <- unname(unlist(strsplit(c,'res\\.')))[2]
    tmp_list[[i]] <- tibble::tibble(pc = pc.use, resolution = r, k_param = k, original_ident_full = list(full_sample_clust[[c]]),
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
