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

sce <- readRDS(file=infile)
expr <- counts(sce)
expr <- matr.filter(as.matrix(expr), min.cells = 10, min.genes = 10)
meta <- as.data.frame(colData(sce))
meta <- meta[, colnames(meta)=='orig.ident' | grepl("^snn", colnames(meta)) | grepl("_res.", colnames(meta))]
meta <- meta[colnames(expr),]
colnames(meta) <- gsub("^snn","walktrap_snn",colnames(meta))
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

plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 60*1024^3)
set.seed(1129)
idx <- which(grepl("snn", colnames(meta)) | grepl("_res.", colnames(meta)))
idx_success <- NULL
rlist <- future_lapply(X = idx, future.seed=TRUE, FUN = function(i) {
        tryCatch({
            res <- rogue(expr, labels = meta[,i], samples = meta$orig.ident, platform = "UMI", span = 0.6)
	    cat(sprintf("%s: ROGUE calculation finished", colnames(meta)[i]),file=logfile,append=T,sep="\n")
	    idx_success <- c(idx_success, i)
	    return(res)
        }, error=function(cond){
	    cat(sprintf("%s: ROGUE calculation failed due to %s", colnames(meta)[i], cond),file=logfile,append=T,sep="\n")
        })
})
rlist <- setNames(rlist, colnames(meta)[idx_success])
save(rlist, file=tmpfile1)

plist <- future_lapply(rlist, future.seed=TRUE, function(r,n){
	r %>% tidyr::gather(key = clusters, value = ROGUE) %>%
	ggplot(aes(clusters, ROGUE)) + geom_boxplot(color = "#FF3E96",outlier.shape = NA) +
	geom_point(color = "#FF3E96", size = 1.5) + theme_bw() + theme(axis.text = element_text(size = 12,
	colour = "black"), axis.title = element_text(size = 13, colour = "black")) + labs(title = n, x = "Clusters", y = "ROGUE")
}, n=names(rlist))
save(plist, file=tmpfile2)

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
