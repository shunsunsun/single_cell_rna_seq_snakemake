library(tidyverse)
library(Seurat)
pc.use=c(10,30,50)
k=20
resolution=c(0.8,1,1.2)
for (method in c("CCA","Hmy","fastMNN")){
    input <- paste0("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_Batch",method,"_clustered.rds")
    se <- readRDS(file=input)
    dat.list <- list()
    i=1
    for(p in pc.use){
        for(r in resolution){
            clustCol <- colnames(se[[]])[grepl(paste0(p,"_snn",k,"_res.",r),colnames(se[[]]))]
            dat.list[[i]] <- tibble::tibble(pc = p, resolution = r, k_param = k, original_ident_full = list(se[[clustCol]]))
            i <- i+1
        }
    }
    gather_idents<- do.call(bind_rows, dat.list)
    saveRDS(gather_idents, file=gsub("_clustered","4clustStab",input))
}
