library(Seurat)

#gej_im <- readRDS(file = "./data/FACS_based/GEJ_IM_pipeCompflted_clustered_SingleRann.rds")
gej_im <- readRDS(file = "./data/FACS_based/ESCC_IM_pipeCompflted_clustered_SingleRann.rds")
#gej_ep <- readRDS(file = "./data/FACS_based/GEJ_EP_pipeCompflted_clustered_SingleRann.rds")
gej_ep <- readRDS(file = "./data/FACS_based/ESCC_EP_pipeCompflted_clustered_SingleRann.rds")

unique(gej_im[["SingleR_cluster_label", drop=T]])
unique(gej_ep[["SingleR_cluster_label", drop=T]])

#---T cells

meta <- gej_im@meta.data
IM_T <- meta[grepl("T cells|Th1 cells", meta$SingleR_cluster_label),]$Cell
IM_nonT <- setdiff(meta$Cell, IM_T)
meta <- gej_ep@meta.data
EP_T <- meta[grepl("CD8", meta$SingleR_cluster_label),]$Cell
EP_nonT <- setdiff(meta$Cell, EP_T)

tcells <- readRDS(file="./data/GEJ_QCed_sctNorm_BatchCCA_clustStab/TCells.rds")
tcells <- rownames(tcells[[]])

length(intersect(tcells, IM_T))
length(intersect(tcells, IM_nonT))
length(intersect(tcells, EP_T))
length(intersect(tcells, EP_nonT))

length(IM_T)
length(IM_nonT)
length(EP_nonT)
length(EP_T)


#--B cells

meta <- gej_im@meta.data
IM_B <- meta[grepl("B-cells|Plasma", meta$SingleR_cluster_label),]$Cell
IM_nonB <- setdiff(meta$Cell, IM_B)
meta <- gej_ep@meta.data
EP_B <- meta[grepl("Plasma", meta$SingleR_cluster_label),]$Cell
EP_nonB <- setdiff(meta$Cell, EP_B)

bcells <- readRDS(file="./data/GEJ_QCed_sctNorm_BatchCCA_clustStab/BCells.rds")
bcells <- rownames(bcells[[]])

length(intersect(bcells, IM_B))
length(intersect(bcells, IM_nonB))
length(intersect(bcells, EP_B))
length(intersect(bcells, EP_nonB))

length(IM_B)
length(IM_nonB)
length(EP_nonB)
length(EP_B)


#--Monocytes

meta <- gej_im@meta.data
IM_M <- meta[grepl("Monocytes|Macrophages", meta$SingleR_cluster_label),]$Cell
IM_nonM <- setdiff(meta$Cell, IM_M)
meta <- gej_ep@meta.data
EP_M <- meta[grepl("Monocytes", meta$SingleR_cluster_label),]$Cell
EP_nonM <- setdiff(meta$Cell, EP_M)

mcells <- readRDS(file="./data/GEJ_QCed_sctNorm_BatchCCA_clustStab/Monocytes.rds")
mcells <- rownames(mcells[[]])

length(intersect(mcells, IM_M))
length(intersect(mcells, IM_nonM))
length(intersect(mcells, EP_M))
length(intersect(mcells, EP_nonM))

length(IM_M)
length(IM_nonM)
length(EP_nonM)
length(EP_M)


#--Epithelials

meta <- gej_im@meta.data
IM_E <- meta[grepl("Epithelial", meta$SingleR_cluster_label),]$Cell
IM_nonE <- setdiff(meta$Cell, IM_E)
meta <- gej_ep@meta.data
EP_E <- meta[grepl("Epithelial", meta$SingleR_cluster_label),]$Cell
EP_nonE <- setdiff(meta$Cell, EP_E)

ecells <- readRDS(file="./data/GEJ_QCed_sctNorm_BatchCCA_clustStab/EpithelialCells.rds")
ecells <- rownames(ecells[[]])

length(intersect(ecells, IM_E))
length(intersect(ecells, IM_nonE))
length(intersect(ecells, EP_E))
length(intersect(ecells, EP_nonE))

length(IM_E)
length(IM_nonE)
length(EP_nonE)
length(EP_E)

#--Stromal

meta <- gej_im@meta.data
IM_S <- meta[grepl("Adipocytes|Neurons", meta$SingleR_cluster_label),]$Cell
IM_nonS <- setdiff(meta$Cell, IM_S)
meta <- gej_ep@meta.data
EP_S <- meta[grepl("Fibroblasts", meta$SingleR_cluster_label),]$Cell
EP_nonS <- setdiff(meta$Cell, EP_S)

scells <- readRDS(file="./data/GEJ_QCed_sctNorm_BatchCCA_clustStab/StromalCells.rds")
scells <- rownames(scells[[]])

length(intersect(scells, IM_S))
length(intersect(scells, IM_nonS))
length(intersect(scells, EP_S))
length(intersect(scells, EP_nonS))

length(IM_S)
length(IM_nonS)
length(EP_nonS)
length(EP_S)
