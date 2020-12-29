load(file="../../data/GEJ_QCed_sctNorm_noBatchCorr_rogueImmun.rda")
res <- do.call("rbind", lapply(1:length(rlist), function(i,all){
	r <- all[[i]]
	r$clust_IQR <- apply(r, 1, function(c){
		return(summary(c,na.rm=T)[5]-summary(c,na.rm=T)[2])
	})
	r$clust_median <- apply(r, 1, median, na.rm=T)
	return(data.frame(clustMethod=names(all)[i],median_IQR=median(r$clust_IQR),min_IQR=min(r$clust_IQR),
		max_IQR=max(r$clust_IQR),median_clust_ROGUE=median(r$clust_median),min_clust_ROGUE=min(r$clust_median),
		max_clust_ROGUE=max(r$clust_median)))
},all=rlist))
res%>%arrange(desc(median_clust_ROGUE),median_IQR)%>%head(n=3)
