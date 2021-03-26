library(Seurat)
library(dplyr)
library(tidyr)


getClustFrac <- function(infile, label){
	se <- readRDS(file=infile)
	meta <- se@meta.data
	colnames(meta)[grepl("_res\\.", colnames(meta))] <- "cluster"
	total <- as.data.frame(table(meta$orig.ident))
	colnames(total)=c("pid","total")
	clust <- as.data.frame(table(meta$orig.ident, meta$cluster))
	colnames(clust)=c("pid","cid","clust_n")
	clust$cid=paste0("Clust",clust$cid)
	clust <- as.data.frame(clust %>% pivot_wider(names_from = cid, values_from = clust_n))
	colnames(clust)=gsub("Clust",paste0(label,"_C"),colnames(clust))
	clust <- merge(clust,total)
	for(i in 2:(ncol(clust)-1)){
		clust[,i]=clust[,i]/clust$total
	}
	return(clust)
}

##ACGEJ
glmpca <- getClustFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustStab/TCells_rc_glmpca10_k150_res0.15.rds", "glmpca")
hmyallhvg <- getClustFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustStab/TCells_sct_allhvg_hmy20_k100_res0.05.rds", "hmyallhvg")
hmy <- getClustFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustStab/TCells_sct_hmy20_k100_res0.05.rds", "hmy")
clinic <- read.table(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/clinicInfo/GEJ_40patients.tsv",header=T)

clust <- merge(merge(glmpca[,-ncol(glmpca)], hmyallhvg[,-ncol(hmyallhvg)]), hmy[,-ncol(hmy)])
clust$site=ifelse(grepl("N",clust$pid),"adjacent","tumor")
frac_cols <- which(grepl("_C", colnames(clust)))
NvsT <- NULL
for(f in frac_cols){
	tmp <- wilcox.test(clust[,f]~clust$site)
	NvsT <- rbind(NvsT, data.frame(clust=colnames(clust)[f], T=median(clust[clust$site=="tumor",f]), 
		N=median(clust[clust$site=="adjacent",f]),p=tmp$p.value))
}
write.table(NvsT[NvsT$p<=0.1, ],col.names=T,row.names=F,sep="\t",quote=F)

clust_T <- clust%>%filter(!grepl("N",pid))
clust_T$pid <- gsub("T","", clust_T$pid)

clinic$Age_group1=ifelse(clinic$Age<=median(clinic$Age), "young", "old")
clinic$Age_group2=ifelse(clinic$Age<median(clinic$Age), "young", "old")
clinic$Stage=ifelse(clinic$Grade<3, "early", "late")
clust_clinic <- merge(clust_T, clinic,by.x="pid",by.y="SampleID")
frac_cols <- which(grepl("_C", colnames(clust_clinic)))
res <- NULL
for(f in frac_cols){
	gender.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Gender)$p.value
	age1.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Age_group1)$p.value
	age2.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Age_group2)$p.value
	smoke.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$smoke)$p.value
	drink.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$drink)$p.value
	stage.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Stage)$p.value
	res <- rbind(res, data.frame(clust=colnames(clust_clinic)[f], gender=gender.p, age1=age1.p, age2=age2.p, smoke=smoke.p, drink=drink.p, stage=stage.p))
}
sig <- NULL
for(i in 1:nrow(res)){
	for(j in 1:ncol(res)){
		if(res[i,j]<=0.05){
			sig=rbind(sig, data.frame(clust=res[i,]$clust,clinic=colnames(res)[j],wilcox_nomial_p=res[i,j]))
		}
	}
}
write.table(sig, col.names=T,row.names=F,sep='\t',quote=F)

#ESCC
glmpca <- getClustFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/TCells_rc_glmpca20_k100_res0.15.rds", "glmpca")
hmyallhvg <- getClustFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/TCells_sct_allhvg_hmy20_k100_res0.05.rds", "hmyallhvg")
hmy <- getClustFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/TCells_sct_hmy20_k100_res0.05.rds", "hmy")
clinic <- read.table(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/clinicInfo/ESCC_60patients.tsv",header=T)

clust <- merge(merge(glmpca[,-ncol(glmpca)], hmyallhvg[,-ncol(hmyallhvg)]), hmy[,-ncol(hmy)])
clust$site=ifelse(grepl("N",clust$pid),"adjacent","tumor")
frac_cols <- which(grepl("_C", colnames(clust)))
NvsT <- NULL
for(f in frac_cols){
	tmp <- wilcox.test(clust[,f]~clust$site)
	NvsT <- rbind(NvsT, data.frame(clust=colnames(clust)[f], T=median(clust[clust$site=="tumor",f]), 
		N=median(clust[clust$site=="adjacent",f]),p=tmp$p.value))
}
write.table(NvsT[NvsT$p<=0.1, ],col.names=T,row.names=F,sep="\t",quote=F)

clust_T <- clust%>%filter(!grepl("N",pid))
clust_T$pid <- gsub("T","", clust_T$pid)

clinic$Age_group1=ifelse(clinic$Age<=median(clinic$Age), "young", "old")
clinic$Age_group2=ifelse(clinic$Age<median(clinic$Age), "young", "old")
clinic$T_stage_group=ifelse(clinic$T_stage %in% c(1,2), "early", "late")
clinic$N_stage_group=ifelse(clinic$N_stage==0, "early", "late")
clinic$P_stage_group1=ifelse(clinic$Pathologic_stage %in% c("I", "II"), "early", "late")
clinic$P_stage_group2=ifelse(clinic$Pathologic_stage=="I", "early", "late")
clust_clinic <- merge(clust_T, clinic,by.x="pid",by.y="SampleID")
frac_cols <- which(grepl("_C", colnames(clust_clinic)))
res <- NULL
for(f in frac_cols){
	gender.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Gender)$p.value
	age1.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Age_group1)$p.value
	age2.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Age_group2)$p.value
	smoke.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Smoking_status)$p.value
	drink.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Drinking_status)$p.value
	tstage.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$T_stage_group)$p.value
	nstage.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$N_stage_group)$p.value
	pstage1.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$P_stage_group1)$p.value
	pstage2.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$P_stage_group2)$p.value
	res <- rbind(res, data.frame(clust=colnames(clust_clinic)[f], gender=gender.p, age1=age1.p, age2=age2.p,
		smoke=smoke.p, drink=drink.p, tstage=tstage.p, nstage=nstage.p, pstage1=pstage1.p, pstage2=pstage2.p))
}

sig <- NULL
for(i in 1:nrow(res)){
	for(j in 1:ncol(res)){
		if(res[i,j]<=0.05){
			sig=rbind(sig, data.frame(clust=res[i,]$clust,clinic=colnames(res)[j],wilcox_nomial_p=res[i,j]))
		}
	}
}
write.table(sig, col.names=T,row.names=F,sep='\t',quote=F)
