library(dplyr)

load(file="xxx_conservedMarkers.rda")
clusts <- names(conserved_markers_list)
for(c in clusts){
       markers <- conserved_markers_list[[c]]
       fold_change <- which(grepl("avg_logFC", colnames(markers)))
       markers$sample_avg_lfc <- apply(markers, 1, function(x){
              return(mean(as.numeric(unname(x[fold_change]))))
       })
       markers %>% select(gene, max_pval, sample_avg_lfc, description) %>% arrange(desc(sample_avg_lfc))
}

load(file="xxx_allMarkers.rda")
top3 <- all_markers %>% filter(p_val_adj<=0.05) %>% group_by(cluster) %>% top_n(3,avg_logFC)


library(Seurat)
library(dplyr)
load(file="./data/GEJ_QCed_sctNorm_BatchCCA_clustStab/Louvain_clust_k100_pc50_res0.4_allMarkers.rda")
top5 <- as.data.frame(all_markers %>% filter(!grepl("^MT-",gene)) %>% filter(!grepl("^RPL",gene)) %>% filter(p_val_adj<=0.05) %>% group_by(cluster) %>% top_n(5,avg_logFC))
top5 <- as.data.frame(all_markers %>% filter(p_val_adj<=0.05) %>% group_by(cluster) %>% top_n(5,avg_logFC))
for(i in unique(top5$cluster)){
	print(paste0("Cluster ", i, ": ", paste(data.frame(top5 %>% filter(cluster==i))$gene,collapse=", ")))
}
