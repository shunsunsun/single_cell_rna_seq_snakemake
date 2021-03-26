indir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/"
se <- readRDS(file=paste0(indir, "TCells_sct_allavg_hmy20_k100_res0.05.rds"))
indir <- "/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustStab/"
se <- readRDS(file=paste0(indir, "TCells_sct_allhvg_hmy20_k100_res0.05.rds"))
meta <- se@meta.data
gr <- c(3,4,5,6,7,8,9)
gr <- c(3,4,5,6)
tcr_group <- meta[meta$SCT_snn_res.0.05 %in% gr,]
table(tcr_group$SCT_snn_res.0.05, tcr_group$orig.ident)
for(i in gr){
    cl <- meta%>%filter(SCT_snn_res.0.05==i)%>%select(orig.ident)
    cl$id <- rownames(cl)
    cl$cid <- apply(cl, 1, function(x){
        unlist(strsplit(x[2],'_'))[2]
    })
    cl=cl[,c(1,3)]
    colnames(cl)=c("pid","cid")
    fwrite(cl,file=paste0(indir, "TCells_sct_allhvg_hmy20_k100_res0.05_c",i,".txt"),col.names=F,row.names=F,sep='\t', quote=F)
}
