##Example: for k in 40 60 80; do for pc in 30 50; do sbatch -p cn_icg -A gaog_g1 --qos=gaogcnicg -N 1 -n 1 -c 20 ./runRscript.sh addClustering.R ./full_sample_preprocess/full_sample_k_${k}_PC_${pc}.rds ${k} ${pc} 20; sleep 1; done; done

suppressPackageStartupMessages({
	library(Seurat)
	library(tidyverse)
	library(future)
})

args  <- commandArgs(trailingOnly=T)
infile=output=args[1]
k <- as.numeric(args[2])
pc.use<- as.numeric(args[3])
resString <- "resolution=c(0.05,0.1)"
reducString <- "reduc=\"harmony\""
nworker <- min(as.numeric(args[4]),length(availableWorkers()))
PreprocessSubsetData_pars<- "variable.features.n = 3000, score.thresh = 1e-5, n.start = 100, nn.eps = 0"
rseed=1129L
seurat_obj<- readRDS(infile)

PreprocessSubsetData<- function(object,
                                variable.features.n = 3000,
                                num.pc = 50,
                                pc.use = NULL,
                                nworker = 8,
                                score.thresh = 1e-5,
                                sig.pc.thresh = 0.05,
                                n.start = 100,
                                nn.eps = 0,
                                resolution = seq(0.6,1.2,by=0.2),
				#resolution = 0.4,
				reduc = "pca",
                                k.param = 30,
				random.seed=1129L,
                                ...){

	
	cat(sprintf("Use %d workders and random seed %d\n",nworker,random.seed))
	plan("multiprocess", workers = nworker)
	options(future.globals.maxSize = 60*1024^3)
	set.seed(random.seed)


	#Enable method = "igraph" to avoid casting large data to a dense matrix
	object <- FindClusters(object, reduction=reduc, n.start=n.start, resolution=resolution,
		random.seed=random.seed, method="igraph",verbose=FALSE)
	
	options(future.globals.maxSize = 500*1024^2) #500M
	plan(sequential)

	return(object)
}


## this is not subsetted data, but the PreprocessSubsetData function can be used as well for any seurat object
seurat_obj<- eval(parse(text=paste("PreprocessSubsetData", "(", "seurat_obj,", "k.param=", k, ",", resString, ",",
	"pc.use=", pc.use, ",", "nworker=", nworker, ",", "random.seed=", rseed, ",", reducString, ",", PreprocessSubsetData_pars, ")")))
saveRDS(seurat_obj,file=output)
