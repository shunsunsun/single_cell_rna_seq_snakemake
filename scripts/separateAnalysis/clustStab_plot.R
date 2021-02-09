##Ref: https://crazyhottommy.github.io/EvaluateSingleCellClustering/mixture_tidy_idents.html

suppressPackageStartupMessages({
	library(scclusteval)
	library(tidyverse)
	library(patchwork)
	library(Seurat)
})

full_rds <- "../../data/ESCC_QCed_sctNorm_BatchHmy_clustStab/gather_full_sample.rds"
sub_rds <- "../../data/ESCC_QCed_sctNorm_BatchHmy_clustStab/gather_subsample.rds"
outplot1 <- "../../plot/clustStab/ESCC_QCed_sctNorm_BatchHmy_scatterPlot1.pdf"
outplot2 <- "../../plot/clustStab/ESCC_QCed_sctNorm_BatchHmy_scatterPlot3.pdf"

subsample_idents <- readRDS(file=sub_rds)
fullsample_idents <- readRDS(file=full_rds)

subsample_idents$original_ident_full=lapply(subsample_idents$original_ident_full, function(x){
	cells <- rownames(x)
	x=as.factor(x[,1])
	names(x) <- cells
	return(x)
})

subsample_idents$recluster_ident=lapply(subsample_idents$recluster_ident, function(x){
	cells <- rownames(x)
	x=as.factor(x[,1])
	names(x) <- cells
	return(x)
})

fullsample_idents$original_ident_full=lapply(fullsample_idents$original_ident_full, function(x){
	cells <- rownames(x)
	x=as.factor(x[,1])
	names(x) <- cells
	return(x)
})

subsample_idents$round <- as.numeric(subsample_idents$round)
subsample_idents=subsample_idents %>% arrange(pc, k_param, resolution,round)

## check for one of the subsample experiment, the cell identities 
## should be the same before and after clustering, but the cluster identities 
## should be different
#subsample_idents %>% mutate(id = row_number()) %>% filter(pc == 30, resolution == 0.8, k_param == 20)
#identical(names(subsample_idents$original_ident_full[[245]]), names(subsample_idents$recluster_ident[[245]]))
#table(subsample_idents$original_ident_full[[245]])
#table(subsample_idents$recluster_ident[[245]])

#Jaccard Raincloud plot for different resolutions
subsample_idents_list<- subsample_idents %>% 
  group_by(pc, resolution, k_param) %>% 
  nest()

#######Skip: Example
#subsample_idents_list %>% ungroup() %>% mutate(id = row_number()) %>%
#  filter(pc == 30, resolution == 0.8, k_param == 20)
#subsample_idents_list$data[[13]]
## for the n times repeating, matching the clusters before and after reclustering
## and assign the jaccard for that cluster
#AssignHighestJaccard(subsample_idents_list$data[[13]]$original_ident_full, subsample_idents_list$data[[13]]$recluster_ident)
#jpeg(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/plot/clustStab/jrdRainPlot_13.jpg", res=300)
#JaccardRainCloudPlot(subsample_idents_list$data[[13]]$original_ident_full,
#                         subsample_idents_list$data[[13]]$recluster_ident) + 
#        geom_hline(yintercept = c(0.6, 0.75), linetype = 2) +
#        xlab("cluster id w/ k=20 res=0.8 pc=30") 
#dev.off()
###########	

#make a Raincloud plot for every combination of the parameters
#subsample_idents_list2 <- subsample_idents_list %>%
#  mutate(plot = map(data, ~JaccardRainCloudPlot(.x$original_ident_full, .x$recluster_ident) + geom_hline(yintercept = c(0.6, 0.75), linetype = 2)))

## to save to disk, give a name for each pdf
#subsample_idents_list2<- mutate(subsample_idents_list2, file_name = 
#	paste0(indir, "plot/clustStab/", prefix, "_PC_",pc, "_", "resolution_",resolution, "_", "k_", k_param, ".pdf"))

# save to disk
#walk2(subsample_idents_list2$file_name, subsample_idents_list2$plot, ggsave, width = 30, height = 20, units="in")


##Skip: Example
#p1<- subsample_idents_list2 %>% 
#  filter(resolution == 0.05, pc == 15, k_param ==20) %>%
#  pull(plot) %>% `[[`(1)  + ggtitle("resolution = 0.05, pc = 15, k.param = 20")


## Skip: example for one set of parameter: k=8, res=0.6, pc = 20
## return a list
#AssignStableCluster(subsample_idents_list$data[[55]]$original_ident,
#                    subsample_idents_list$data[[55]]$recluster_ident,
#                    jaccard_cutoff = 0.8,
#                    method = "jaccard_percent", 
#                   percent_cutoff = 0.8)
					
#Assign stable clusters					
stable_clusters<- subsample_idents_list %>%
  mutate(stable_cluster = map(data, ~ AssignStableCluster(.x$original_ident_full,
  .x$recluster_ident, jaccard_cutoff = 0.8, method = "jaccard_percent", percent_cutoff = 0.8)))					

#plot scatter plot for different parameters sets				
p1 <- ParameterSetScatterPlot(stable_clusters = stable_clusters,
                        fullsample_idents = fullsample_idents,
                        x_var = "k_param",
                        y_var = "number",
                        facet_rows = "resolution",
                        facet_cols = "pc")
ggsave(p1, file= outplot1, width=4, height=5, units="in")			
					
#p2=ParameterSetScatterPlot(stable_clusters = stable_clusters,
#                        fullsample_idents = fullsample_idents,
#                        x_var = "resolution",
#                        y_var = "number",
#                        facet_rows = "k_param",
#                        facet_cols = "pc")
#ggsave(p2, file= outplot2, width=4, height=5, units="in")						
					
#Calculate percentage of cells in stable clusters
stable_and_full<- left_join(stable_clusters, fullsample_idents)
stable_and_full<- stable_and_full %>% 
  mutate(precentage_in_stable = map2_dbl(original_ident_full, stable_cluster, 
  function(x, y) CalculatePercentCellInStable(x,y$stable_cluster)))
  
p3 <- ParameterSetScatterPlot(stable_clusters = stable_clusters,
                        fullsample_idents = fullsample_idents,
                        x_var = "k_param",
                        y_var = "percentage",
                        facet_rows = "resolution",
                        facet_cols = "pc") + ggtitle("percentage of cells in stable clusters")
ggsave(p3, file= outplot2, width=4, height=5, units="in")
