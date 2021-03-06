suppressPackageStartupMessages({
	library(Seurat)
	library(tidyverse)
	library(future)
})

## see https://bitbucket.org/snakemake/snakemake/issues/917/enable-stdout-and-stderr-redirection
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

seurat_obj<- readRDS(snakemake@input[[1]])
k <- snakemake@wildcards[["k"]]
pc.use<- snakemake@wildcards[["pc"]]
resString <- as.character(snakemake@params[["res"]])
resVal <- as.numeric(unname(unlist(strsplit(resString,','))))
if(length(resVal)>1){
	resString=paste0("c(",resString,")")
}
print(paste0(resString, ":", resVal))
reducString <- snakemake@params[["reduc"]]
nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
output <- snakemake@output[[1]]
PreprocessSubsetData_pars<- snakemake@params[["PreprocessSubsetData_pars"]]
rseed=1129L

PreprocessSubsetData<- function(object,
                                variable.features.n = 3000,
                                num.pc = 50,
                                pc.use = NULL,
                                nworker = 8,
                                score.thresh = 1e-5,
                                sig.pc.thresh = 0.05,
                                n.start = 100,
                                nn.eps = 0,
                                resolution = seq(0.4,0.8,by=0.2),
				#resolution = 0.4,
				reduc = "pca",
                                k.param = 30,
				random.seed=1129L,
                                ...){

	
	cat(sprintf("Use %d workders and random seed %d\n",nworker,random.seed))
	plan("multiprocess", workers = nworker)
	options(future.globals.maxSize = 60*1024^3)
	set.seed(random.seed)


	#https://satijalab.org/seurat/v3.0/integration.html
	if(reduc == "pca" && !("pca" %in% names(object@reductions))){
		object <- RunPCA(object, verbose = FALSE, npcs=num.pc, features = VariableFeatures(object = object)) ##the default value of npcs = the default value of num.pc = 50
	}
	if(reduc == "glmpca"){
		reduc <- paste0(reduc,pc.use)
		print(paste0("reduc=",reduc))
	}

#       if (is.null(pc.use)){
#               object<- JackStraw( object = object, num.replicate = 100, dims = num.pc)
#               object <- ScoreJackStraw(object = object, dims = 1:num.pc, score.thresh = score.thresh)
#               PC_pvalues<- object@reductions$pca@jackstraw@overall.p.values
#               ## determin how many PCs to use
#               pc.use<- min(which(PC_pvalues[,"Score"] > sig.pc.thresh)) -1
#       }
	
	object <- FindNeighbors(object,reduction=reduc,dims=1:pc.use,k.param=k.param,verbose=FALSE,force.recalc=TRUE)
	
	#Enable method = "igraph" to avoid casting large data to a dense matrix
	object <- FindClusters(object, n.start=n.start, resolution=resolution, random.seed=random.seed, method="igraph",verbose=FALSE)
	
	options(future.globals.maxSize = 500*1024^2) #500M
	plan(sequential)

	return(object)
}


## this is not subsetted data, but the PreprocessSubsetData function can be used as well for any seurat object
command <- paste("PreprocessSubsetData", "(", "seurat_obj,", "k.param=", k, ",", "resolution=", resString, ",",
        "pc.use=", pc.use, ",", "nworker=", nworker, ",", "random.seed=", rseed, ",", reducString, ",", PreprocessSubsetData_pars, ")")
print(paste0("Command=",command))
seurat_obj<- eval(parse(text=command))
saveRDS(seurat_obj,file=output)
