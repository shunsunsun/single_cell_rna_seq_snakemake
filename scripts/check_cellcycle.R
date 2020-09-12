input <- snakemake@input[[1]] ## seurat object as stored as *.rds
filename <- unname(unlist(strsplit(basename(input),"_")))
outfile=snakemake@output[["split"]]
outdir <- dirname(outfile)
if(!dir.exists(outdir)){
        dir.create(outdir,recursive=T)
}
outfile2=snakemake@output[["overlay"]]

library(Seurat)

## Read input
seurat_obj <- readRDS(input)

# Normalize the counts
if(snakemake@params[["sctPreNorm"]]){
	seurat_phase <- SCTransform(seurat_obj, vars.to.regress = "subsets_Mito_percent", verbose = FALSE)
}else{
	seurat_phase <- NormalizeData(seurat_obj)
}


# Load cell cycle markers
load(snakemake@params[["ccgenes"]])

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, g2m.features = g2m_genes, s.features = s_genes)


# (cannot use in our HPC) View cell cycle scores and phases assigned to cells                                 
# View(seurat_phase@meta.data)


# Identify the most variable genes (default setting)
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)


# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
jpeg(outfile, width = 18, height = 5, units="in", res=300)
pca=DimPlot(seurat_phase, reduction = "pca", group.by= "Phase", split.by = "Phase")
print(pca)
dev.off()

jpeg(outfile2, width = 8, height = 6, units="in", res=300)
pca=DimPlot(seurat_phase, reduction = "pca", group.by= "Phase")
print(pca)
dev.off()




