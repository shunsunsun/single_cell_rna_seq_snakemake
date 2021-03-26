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
getThFrac <- function(infile, singleRannt, cl){
	se <- readRDS(file=infile)
	load(file=singleRannt)
	se[["SingleR"]]=cell.pred$labels
	meta <- se@meta.data
	colnames(meta)[grepl("_res\\.", colnames(meta))] <- "cluster"
	total <- as.data.frame(table(meta$orig.ident))
	colnames(total)=c("pid","total")
	meta <- meta[meta$cluster==cl & (grepl("CD4", meta$SingleR)|grepl("Th", meta$SingleR)|grepl("helper", meta$SingleR)),]
	clust <- as.data.frame(table(meta$orig.ident))
	colnames(clust)=c("pid","clust_n")
	clust <- merge(clust,total)
	clust$frac=clust$clust_n/clust$total
	return(clust)
}

##ACGEJ
glmpca <- getClustFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustStab/TCells_rc_glmpca10_k150_res0.15.rds", "glmpca")
hmyallhvg <- getClustFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustStab/TCells_sct_allhvg_hmy20_k100_res0.05.rds", "hmyallhvg")
hmy <- getClustFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/GEJ_QCed_sctNorm_BatchCCA_clustStab/TCells_sct_hmy20_k100_res0.05.rds", "hmy")
clinic <- read.table(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/clinicInfo/GEJ_40patients.tsv",header=T)

clust <- merge(merge(glmpca[,-ncol(glmpca)], hmyallhvg[,-ncol(hmyallhvg)]), hmy[,-ncol(hmy)])
#clust <- glmpca[,-ncol(glmpca)]
clust$site=ifelse(grepl("N",clust$pid),"adjacent","tumor")
frac_cols <- which(grepl("_C", colnames(clust)))
#frac_cols <- which(grepl("frac", colnames(clust)))
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
#frac_cols <- which(grepl("frac", colnames(clust)))
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
p_thresh=0.05
sig <- NULL
for(i in 1:nrow(res)){
	for(j in 1:ncol(res)){
		if(res[i,j]<=p_thresh){
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
#clust <- glmpca[,-ncol(glmpca)]
clust <- getThFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/TCells_rc_glmpca20_k100_res0.15.rds",
	"/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/ESCC_QCed_sctNorm_BatchHmy_clustStab/TCells_rc_glmpca20_k100_res0.15_singleR_perCell_mi.rda",3)
clust$site=ifelse(grepl("N",clust$pid),"adjacent","tumor")
frac_cols <- which(grepl("_C", colnames(clust)))
#frac_cols <- which(grepl("frac", colnames(clust)))
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
clust_clinic <- merge(clust_T, clinic, by.x="pid",by.y="SampleID")
frac_cols <- which(grepl("_C", colnames(clust_clinic)))
#frac_cols <- which(grepl("frac", colnames(clust_clinic)))
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

p_thresh=0.1
sig <- NULL
for(i in 1:nrow(res)){
	for(j in 1:ncol(res)){
		if(res[i,j]<=p_thresh){
			sig=rbind(sig, data.frame(clust=res[i,]$clust,clinic=colnames(res)[j],wilcox_nomial_p=res[i,j]))
		}
	}
}
write.table(sig, col.names=T,row.names=F,sep='\t',quote=F)

epi_program <- readRDS(file="./data/SampleInfo_EpiProgram_0.7.rds")
clust_epi <- merge(clust_T, epi_program, by.x="pid",by.y="sample")
frac_cols <- which(grepl("_C", colnames(clust_epi)))
prog_cols <- (ncol(clust_epi)-7):ncol(clust_epi)
res <- NULL
for(f in frac_cols){
    for (p in prog_cols){
        tmp <- cor.test(clust_epi[,f], clust_epi[,p], method="spearman")
        res <- rbind(res, data.frame(clust=colnames(clust_epi)[f], program=colnames(clust_epi)[p], p=tmp$p.value, rho=tmp$estimate))
    }
}
p_thresh=0.05
write.table(res[res$p<=p_thresh, ], col.names=T,row.names=F,sep='\t',quote=F)

p = ggscatter(clust_epi,x = "glmpca_C3",y = "Epi1",size = 1,ncol = 5,xlab = 'TCells_C3',ylab = 'Epi1 Score',conf.int = TRUE) +
    geom_smooth(method = "lm",se = FALSE, show.legend = FALSE, color = "blue",size = 0.5) +
    stat_cor(method = "spearman", label.x=0.4) + theme(aspect.ratio = 1) +
    theme(panel.background  = element_rect(fill = 'white'),axis.ticks = element_line(colour =  'black'),
      axis.line = element_line(colour =  'black'),axis.text = element_text(colour = 'black'),legend.position = "none")
ggsave(plot=p, file="plot/oralReport_0322/ESCC_TCellsC3_Epi1Score.pdf", width = 5, height = 5, units="in",dpi=300)

p = ggscatter(clust_epi,x = "glmpca_C2",y = "Epi1",size = 1,ncol = 5,xlab = 'BCells_C2',ylab = 'Epi1 Score',conf.int = TRUE) +
    geom_smooth(method = "lm",se = FALSE, show.legend = FALSE, color = "blue",size = 0.5) +
    stat_cor(method = "spearman", label.x=0.2) + theme(aspect.ratio = 1) +
    theme(panel.background  = element_rect(fill = 'white'),axis.ticks = element_line(colour =  'black'),
      axis.line = element_line(colour =  'black'),axis.text = element_text(colour = 'black'),legend.position = "none")
ggsave(plot=p, file="plot/oralReport_0322/ESCC_BCellsC2_Epi1Score.pdf", width = 5, height = 5, units="in",dpi=300)


p1 <- ggboxplot(clust_epi,x="T_stage_group",y="Epi1", color="T_stage_group", palette = "jco", add = "jitter", ylab="Epi1 Score", size=1.2) + stat_compare_means(method = "wilcox.test",size=7, label.x=1.3, label.y=1) + font("x.text", size = 24, face="bold")
p2 <- ggboxplot(clust_epi,x="P_stage_group2",y="Epi1",color="P_stage_group2", palette = "jco", add = "jitter", ylab="Epi1 Score", size=1.2) + stat_compare_means(method = "wilcox.test",size=7,label.x=1.3, label.y=1) + font("x.text", size = 24, face="bold")
p=p1+p2
ggsave(plot=p,file="plot/oralReport_0322/ESCC_Epi1_stage.pdf",width = 12, height = 8, units="in",dpi=300)

##Combined
glmpca <- getClustFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/Combined/TCells_rc_glmpca10_k100_res0.15.rds", "glmpca")
hmyallhvg <- getClustFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/Combined/TCells_sct_allhvg_hmy20_k100_res0.05.rds", "hmyallhvg")
hmy <- getClustFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/Combined/TCells_sct_hmy20_k100_res0.05.rds", "hmy")
clinic_g <- read.table(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/clinicInfo/GEJ_40patients.tsv",header=T)
clinic_g$cancertype=rep("G",n=nrow(clinic_g))
clinic_e <- read.table(file="/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/clinicInfo/ESCC_60patients.tsv",header=T)
clinic_e$cancertype=rep("E",n=nrow(clinic_e))
clinic_g$T_stage <- ifelse(grepl("T4", clinic_g$TNM), 4, 
    ifelse(grepl("T3", clinic_g$TNM), 3, 
    ifelse(grepl("T2", clinic_g$TNM), 2, 
    ifelse(grepl("T1", clinic_g$TNM), 1, NA))))
clinic_g$N_stage <- ifelse(grepl("N4", clinic_g$TNM), 4, 
    ifelse(grepl("N3", clinic_g$TNM), 3, 
    ifelse(grepl("N2", clinic_g$TNM), 2, 
    ifelse(grepl("N1", clinic_g$TNM), 1, 
    ifelse(grepl("N0", clinic_g$TNM), 0, NA)))))
colnames(clinic_g)=gsub("smoke","Smoking_status",colnames(clinic_g))
colnames(clinic_g)=gsub("drink","Drinking_status",colnames(clinic_g))
clinic_g$Pathologic_stage=ifelse(clinic_g$Grade==1, "I", ifelse(clinic_g$Grade==2, "II", ifelse(clinic_g$Grade==3, "III", "IV")))
clinic_g$Gender=ifelse(clinic_g$Gender==1, "Male", "Female")
commonCols <- intersect(colnames(clinic_e), colnames(clinic_g))
clinic_g <- clinic_g[,commonCols]
clinic_e <- clinic_e[,commonCols]
clinic <- rbind(clinic_g, clinic_e)

clinic$Age_group1=ifelse(clinic$Age<=median(clinic$Age), "young", "old")
clinic$Age_group2=ifelse(clinic$Age<median(clinic$Age), "young", "old")
clinic$T_stage_group1=ifelse(clinic$T_stage %in% c(1,2), "early", "late")
clinic$T_stage_group2=ifelse(clinic$T_stage==1, "early", "late")
clinic$N_stage_group1=ifelse(clinic$N_stage==0, "early", "late")
clinic$N_stage_group2=ifelse(clinic$N_stage %in% c(1,2), "early", "late")
clinic$P_stage_group1=ifelse(clinic$Pathologic_stage %in% c("I", "II"), "early", "late")
clinic$P_stage_group2=ifelse(clinic$Pathologic_stage=="I", "early", "late")

###Check potential bias in the clinical features of ACGEJ and ESCC
table(clinic$cancertype, clinic$T_stage_group1)
fisher.test(table(clinic$cancertype, clinic$P_stage_group2))
fisher.test(table(clinic$cancertype, clinic$P_stage_group2),alternative="greater")

clust <- merge(merge(glmpca[,-ncol(glmpca)], hmyallhvg[,-ncol(hmyallhvg)]), hmy[,-ncol(hmy)])
#clust <- glmpca[,-ncol(glmpca)]
clust <- getThFrac("/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/Combined/TCells_rc_glmpca10_k100_res0.15.rds","/gpfs2/gaog_pkuhpc/users/liny/GEJ_singleCell/data/Combined/TCells_rc_glmpca10_k100_res0.15_singleR_perCell_mi.rda",3)

frac_cols <- which(grepl("_C", colnames(clust)))
#frac_cols <- which(grepl("frac", colnames(clust)))
clust$site=ifelse(grepl("N",clust$pid),"adjacent","tumor")
NvsT <- NULL
for(f in frac_cols){
	tmp <- wilcox.test(clust[,f]~clust$site)
	NvsT <- rbind(NvsT, data.frame(clust=colnames(clust)[f], T=median(clust[clust$site=="tumor",f]), 
		N=median(clust[clust$site=="adjacent",f]),p=tmp$p.value))
}
write.table(NvsT[NvsT$p<=0.1, ],col.names=T,row.names=F,sep="\t",quote=F)

clust_T <- clust%>%filter(!grepl("N",pid))
clust_T$pid <- gsub("T","", clust_T$pid)

clust_clinic <- merge(clust_T, clinic,by.x="pid",by.y="SampleID")
frac_cols <- which(grepl("_C", colnames(clust_clinic)))
#frac_cols <- which(grepl("frac", colnames(clust)))
res <- NULL
for(f in frac_cols){
	gender.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Gender)$p.value
	age1.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Age_group1)$p.value
	age2.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Age_group2)$p.value
	smoke.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Smoking_status)$p.value
	drink.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$Drinking_status)$p.value
	tstage1.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$T_stage_group1)$p.value
    tstage2.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$T_stage_group2)$p.value
	nstage1.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$N_stage_group1)$p.value    
	nstage2.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$N_stage_group2)$p.value
	pstage1.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$P_stage_group1)$p.value
	pstage2.p <- wilcox.test(clust_clinic[,f] ~ clust_clinic$P_stage_group2)$p.value
	res <- rbind(res, data.frame(clust=colnames(clust_clinic)[f], gender=gender.p, age1=age1.p, age2=age2.p,
		smoke=smoke.p, drink=drink.p, tstage1=tstage1.p, tstage2=tstage2.p, nstage1=nstage1.p, 
        nstage2=nstage2.p, pstage1=pstage1.p, pstage2=pstage2.p))
}

p_thresh=0.05
sig <- NULL
for(i in 1:nrow(res)){
	for(j in 1:ncol(res)){
		if(res[i,j]<=p_thresh){
			sig=rbind(sig, data.frame(clust=res[i,]$clust,clinic=colnames(res)[j],wilcox_nomial_p=res[i,j]))
		}
	}
}
write.table(sig, col.names=T,row.names=F,sep='\t',quote=F)

median(clust_clinic[clust_clinic$P_stage_group1=="early", ]$glmpca_C0)
median(clust_clinic[clust_clinic$P_stage_group1=="late", ]$glmpca_C0)


