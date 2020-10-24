suppressPackageStartupMessages(library(Seurat))

vstpca <- snakemake@input[[1]]
pcaDrName <- snakemake@params[["pca_dr"]]
glmPrefix <- snakemake@params[["glmpca"]]
glmDrName <- snakemake@params[["glm_dr"]]
outfile <- snakemake@output[[1]]

se1 <- readRDS(file=vstpca)
se <- DietSeurat(se1, counts = TRUE, data = TRUE, scale.data = TRUE, dimreducs = NULL)
se[["pca"]]=se1[[pcaDrName]]
Key(se[["pca"]])="PC_"
for(d in c(10,20)){
	se2 <- readRDS(file=paste0(glmPrefix,"_ndim",d,".rds"))
	se[[paste0("glmpca",d)]]=se2[[glmDrName]]
	Key(se[[paste0("glmpca",d)]])=paste0("n",d,"GLMPC_")
}
saveRDS(se,file=outfile)
#if(file.exists(outfile)){
#	file.remove(vstpca)
#	file.remove(paste0(glmPrefix,"_ndim10.rds"))
#	file.remove(paste0(glmPrefix,"_ndim20.rds"))
#}
