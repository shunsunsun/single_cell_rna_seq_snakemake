library("edgeR",lib="/gpfs2/gaog_pkuhpc/users/liny/tools/R-3.6.0-library/")
lost <- calculateAverage(counts(sce)[,sce$discard])
kept <- calculateAverage(counts(sce)[,!sce$discard])
logged <- cpm(cbind(lost, kept), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2] 
abundance <- rowMeans(logged)
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16) 
points(abundance[is.mito], logFC[is.mito], col="dodgerblue", pch=16) 

library(org.Hs.eg.db,lib="/gpfs2/gaog_pkuhpc/users/liny/tools/R-3.6.0-library/")
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
platelet <- c("PF4", "PPBP", "CAVIN2")
ids <- mapIds(org.Hs.eg.db, keys=platelet, column="ENSEMBL", keytype="SYMBOL")
points(abundance[ids], logFC[ids], col="orange", pch=16)
