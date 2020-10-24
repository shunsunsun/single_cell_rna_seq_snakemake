library(Seurat)
library(ggplot2)
library(patchwork)
# set up future for parallelization
library(future)
library(future.apply)

input <- snakemake@input[[1]] ## *_flted.rds
out_seurat <- snakemake@output[["rds"]]
pca_plot=snakemake@output[["pca"]]
tsne_plot=snakemake@output[["tsne"]]
umap_plot=snakemake@output[["umap"]]
combine_plot=snakemake@output[["combine"]]
nPC <- as.numeric(snakemake@params[["nPC"]])
nworker <- min(as.numeric(snakemake@threads),length(availableWorkers()))
cat(sprintf("Use %d workders\n",nworker))

outdir1 <- dirname(out_seurat)
if(!dir.exists(outdir1)){
        dir.create(outdir1,recursive=T)
}
outdir2 <- dirname(pca_plot)
if(!dir.exists(outdir2)){
        dir.create(outdir2,recursive=T)
}
#if(snakemake@params[["regressNum"]]){
#    regVars=c('subsets_Mito_percent', 'nFeature_RNA', 'nCount_RNA')
#else{
    regVars=c('subsets_Mito_percent')
#}


## Read input
flted_seurat <- readRDS(input)

## Load cell cycle markers
load(snakemake@params[["ccgenes"]])


plan("multiprocess", workers = nworker)

## Conduct normalization
if(!snakemake@params[["normPerSample"]]){
	if(!snakemake@params[["sctPreNorm"]]){
        	seurat_phase <- NormalizeData(flted_seurat)
        	seurat_phase <- CellCycleScoring(seurat_phase, g2m.features=g2m_genes, s.features=s_genes)
        	norm_seurat <- SCTransform(seurat_phase, vars.to.regress = c(regVars, 'S.Score', 'G2M.Score'))
	}else{
        	seurat_phase <- SCTransform(flted_seurat, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = regVars)
        	seurat_phase <- CellCycleScoring(seurat_phase, s.features = s_genes, g2m.features = g2m_genes, assay = 'SCT', set.ident = TRUE)
        	norm_seurat <- SCTransform(seurat_phase, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c(regVars, 'S.Score', 'G2M.Score'))
	}

}else{
	options(future.globals.maxSize = 20*1024^3) #20G
	if(!file.exists(gsub("rds","rda",out_seurat))){
	
		# Split seurat object to perform cell cycle scoring and SCTransform on all samples
		split_seurat <- SplitObject(flted_seurat, split.by = "orig.ident")

		# Normalize each sample
		split_seurat <- future_lapply(X = split_seurat, FUN = function(s) {
			print(unique(s@meta.data$orig.ident))
			##After the filtering, sample-level cell/feature number criterion may not be met
			x=CreateSeuratObject(counts=GetAssayData(object=s, slot='counts'), meta.data=s@meta.data, min.cells=3, min.features=200)
			if(!snakemake@params[["sctPreNorm"]]){
	        		x <- tryCatch({
					tmp <- NormalizeData(x)
					tmp <- CellCycleScoring(tmp, g2m.features=g2m_genes, s.features=s_genes)
	        			return(SCTransform(tmp, vars.to.regress = c(regVars, 'S.Score', 'G2M.Score')))
				}, error=function(cond){
					#if(grepl("Insufficient data values to produce", cond, fixed=TRUE)){
					return(SCTransform(x, vars.to.regress = regVars))
					#}
				})
			}else{
				x <- tryCatch({
					tmp <- SCTransform(x, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = regVars)
	        			tmp <- CellCycleScoring(tmp, s.features = s_genes, g2m.features = g2m_genes, assay = 'SCT', set.ident = TRUE)
	        			return(SCTransform(tmp, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c(regVars, 'S.Score', 'G2M.Score')))
				}, error=function(cond){
					return(SCTransform(x, vars.to.regress = regVars))
				})
			}
		})

		# Select the most variable features to use for integration
		integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)

		# Prepare the SCT list object for integration
		split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features)

		save(split_seurat, integ_features, file=gsub("rds","rda",out_seurat))

	}else{
		load(file=gsub("rds","rda",out_seurat))
	}
	# Find best buddies - can take a while to run
	if(snakemake@params[["reduction"]]=="cca"){
		integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = "SCT", anchor.features = integ_features, verbose = TRUE, dims=1:nPC)
	}else{
		split_seurat <- future_lapply(X = split_seurat, FUN = RunPCA, verbose=FALSE, features=integ_features)
		integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = "SCT", anchor.features = integ_features, reduction = "rpca", verbose = TRUE, dims=1:nPC)
	}
	# Integrate across conditions	
	norm_seurat <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", dims=1:nPC)
	options(future.globals.maxSize = 500*1024^2) #500M	
}

saveRDS(norm_seurat, out_seurat)

# Perform dimensionality reduction by PCA, TSNE and UMAP embedding
set.seed(1129)
norm_seurat <- FindVariableFeatures(norm_seurat,selection.method = "vst", nfeatures = 2000)
norm_seurat <- RunPCA(norm_seurat,features=VariableFeatures(norm_seurat))
saveRDS(norm_seurat, out_seurat)

norm_seurat <- RunUMAP(norm_seurat, dims = 1:nPC)
norm_seurat <- RunTSNE(norm_seurat, dims = 1:nPC)

norm_seurat <- FindNeighbors(norm_seurat, dims = 1:nPC)
norm_seurat <- FindClusters(norm_seurat)

plan(sequential)

# Visualization

# group_by_sample
plot1 <- DimPlot(norm_seurat, reduction = "pca", group.by="orig.ident") + NoLegend()
# group_by_cluster
plot2 <- DimPlot(norm_seurat, reduction = "pca", label=T) + NoLegend()
plot3 <- ElbowPlot(norm_seurat, ndims=30, reduction="pca") 
p = plot1 + plot2 + plot3
ggsave(file=pca_plot, plot = p, width = 18, height = 5, device="jpeg",units="in",dpi=300)

#group_by_cluster
plot1 = DimPlot(norm_seurat, reduction = "tsne", label=T) + NoLegend()
#group_by_sample
plot2 = DimPlot(norm_seurat, reduction = "tsne", group.by='orig.ident') + NoLegend() 
g <- plot1 + plot2
ggsave(file=tsne_plot, plot = g, width = 10, height = 5)

#group_by_cluster
plot3 = DimPlot(norm_seurat, reduction = "umap", label=T) + NoLegend()
#group_by_sample
plot4 = DimPlot(norm_seurat, reduction = "umap", group.by='orig.ident') + NoLegend()
g <- plot3+plot4
ggsave(file=umap_plot, plot = g, width = 10, height = 5)

#put tsne and umap side by side
g <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave(file=combine_plot, plot = g, width = 10, height = 5)
