#devtools::install_github("PaulingLiu/ROGUE")
suppressMessages({
	library(Seurat)
	library(future.apply)
	library(scater)
	library(tidyverse)
	library(ggplot2)
	library(gridExtra)
	library(ROGUE)
})
args  <- commandArgs(trailingOnly=T)
infile=args[1]
#infile <- snakemake@input[[1]] #GEJ_QCed_*_walktrap_clustered.rds
tmpfile1 <- gsub("walktrap_clustered.rds","rogue.rda",infile)
tmpfile2 <- gsub("walktrap_clustered.rds","roguePlots.rda",infile)
logfile <- paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/log/",gsub("walktrap_clustered.rds","rogue.log",basename(infile)))
nworker=min(as.numeric(args[2]),length(availableWorkers()))
cat(sprintf("Use %d workers", nworker),file=logfile,append=T,sep="\n")
plotDir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/batchCorrect/"
useImmun=args[3]
if(useImmun=='T'){
	tmpfile1 <- gsub("rogue","rogueImmun",tmpfile1)
	tmpfile2 <- gsub("rogue","rogueImmun",tmpfile2)
	logfile <- gsub("rogue","rogueImmun",logfile)
}

sce <- readRDS(file=infile)
if(useImmun=='T'){
    meta <- colData(sce)
    sce <- sce[,grepl("-I",rownames(meta))]
}
if(!grepl("BatchCCA", infile)){
	expr <- counts(sce)
}else{
	expr <- logcounts(sce)
}

plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
set.seed(1129)

if(!file.exists(tmpfile1)){
	expr <- matr.filter(as.matrix(expr), min.cells = 10, min.genes = 10)
	meta <- as.data.frame(colData(sce))
	meta <- meta[, colnames(meta)=='orig.ident' | grepl("snn", colnames(meta)) | grepl("_res.", colnames(meta))]
	meta <- meta[colnames(expr),]
	walktrap_res <- which(grepl("snn",colnames(meta)) & !grepl("_res",colnames(meta)))
	colnames(meta)[walktrap_res]=gsub("snn","walktrap_snn",colnames(meta)[walktrap_res])
	meta$orig.ident=sapply(strsplit(meta$orig.ident,"-"),"[",1)

	#seplotfile <- paste0(plotDir,gsub(".rda","_SEplot.jpg",basename(outfile)))
	#if(!file.exists(seplotfile)){
	#	ent.res <- SE_fun(expr)
	#	jpeg(file=paste0(plotDir,gsub(".rda","_SEplot.jpg",basename(outfile))), width = 8, height = 6, units="cm", res=300)
	#	SEplot(ent.res)
	#	dev.off()
		#rogue.value <- CalculateRogue(ent.res, platform = "UMI")
		#print(paste0("Global ROGUE value: ", rogue.value))
	#}

	idx <- which(grepl("snn", colnames(meta)) | grepl("_res.", colnames(meta)))
	rlist <- future_lapply(idx, future.seed=TRUE, FUN = function(i) {
        	tryCatch({
	            res <- rogue(expr, labels = meta[,i], samples = meta$orig.ident, platform = "UMI", span = 0.6)
		    colnames(res) <- gsub("^","cluster",colnames(res))
		    res <- data.frame(clustMethod=colnames(meta)[i],res)
		    cat(sprintf("%s: ROGUE calculation finished", colnames(meta)[i]),file=logfile,append=T,sep="\n")
		    return(res)
        	}, error=function(cond){
		    cat(sprintf("%s: ROGUE calculation failed due to %s", colnames(meta)[i], cond),file=logfile,append=T,sep="\n")
		    return(NULL)
	        })
	})
	rlist <- Filter(Negate(is.null), rlist)
	for (k in 1:length(rlist)){
		r <- rlist[[k]]
		names(rlist)[k] <- unique(r$clustMethod)
		rlist[[k]] <- r[,colnames(r)!="clustMethod"]	
	}
	save(rlist, file=tmpfile1)
}else{
	load(file=tmpfile1)
	print(paste0("loading ", tmpfile1))
}

if(!file.exists(tmpfile2)){
	plist <- future_lapply(seq_along(rlist), future.seed=TRUE, function(r,n,i){
		r[[i]] %>% tidyr::gather(key = clusters, value = ROGUE) %>% mutate(clusters = gsub("cluster", "", clusters) %>%
		ggplot(aes(clusters, ROGUE)) + geom_boxplot(color = "#FF3E96",outlier.shape = NA) +
		geom_point(color = "#FF3E96", size = 1.5) + theme_bw() + theme(axis.text = element_text(size = 12,
		colour = "black"), axis.title = element_text(size = 13, colour = "black")) + labs(title = n[[i]], x = "Clusters", y = "ROGUE")
	}, r=rlist, n=names(rlist))
	save(plist, file=tmpfile2)
}else{
	load(file=tmpfile2)
        print(paste0("loading ", tmpfile2))
}

plotfile <- paste0(plotDir, gsub("rda","jpg",basename(tmpfile1)))
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

options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
