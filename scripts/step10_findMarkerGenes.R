suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(future))

infile <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]
outdir <- snakemake@params[[1]]
if(!dir.exists(outdir)){
        dir.create(outdir,recursive=T)
}

se <- readRDS(file=infile)
nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))
plan("multiprocess", workers = nworker)
options(future.globals.maxSize = 10*1024^3) #10G
randSeed <- 1129L
set.seed(randSeed)

outPrefix=paste0(outdir,"/",gsub("clustered.rds","",basename(infile)))
for(i in unique(se@meta.data$seurat_clusters)){
  markers_df <- FindMarkers(object = se, ident.1 = i, min.pct = 0.25)
  print(x = head(markers_df))
  markers_genes = rownames(head(x = markers_df, n = 5))
  VlnPlot(object = se,features =markers_genes,log =T )
  ggsave(filename=paste0(outPrefix,'VlnPlot_Cluster',i,'_markers_heatmap.pdf'))
  FeaturePlot(object = se, features=markers_genes )
  ggsave(filename=paste0(outPrefix,'_FeaturePlot_Cluster',i,'_markers_heatmap.pdf'))
}

markers <- FindAllMarkers(object = se, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
save(markers, outfile)
options(future.globals.maxSize = 500*1024^2) #500M
plan(sequential)
