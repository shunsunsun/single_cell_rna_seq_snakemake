  library(ggrepel)
  require(qdapTools)
  require(REdaS)
  # 
  metadata <- phe
  # Keep only cells from tissues that are not brain or pleura 
  metadata <- metadata[-which(metadata$biopsy_site=="Brain" | metadata$biopsy_site=="Pleura"),]
  # Convert to factor with logical order 
  metadata$analysis <- factor(metadata$analysis, levels = c("naive", "grouped_pr", "grouped_pd"))
  # Create table and keep selected cell types 
  meta.temp <- metadata[,c("immuSub", "analysis")]
  # Loop over treatment response categories 
  # Create list to store frequency tables 
  prop.table.error <- list()
  for(i in 1:length(unique(meta.temp$analysis))){
    vec.temp <- meta.temp[meta.temp$analysis==unique(meta.temp$analysis)[i],"immuSub"]
    # Convert to counts and calculate 95% CI 
    # Store in list 
    table.temp <- freqCI(vec.temp, level = c(.95))
    prop.table.error[[i]] <- print(table.temp, percent = TRUE, digits = 3)
    # 
  }
  # Name list 
  names(prop.table.error) <- unique(meta.temp$analysis)

 # Convert to data frame 
  tab.1 <- as.data.frame.array(do.call(rbind, prop.table.error))
  # Add analysis column 
  b <- c()
  a <- c()
  for(i in names(prop.table.error)){
    a <- rep(i,nrow(prop.table.error[[1]]))
    b <- c(b,a)
  }
  tab.1$analysis <- b
  # Add common cell names 
  aa <- gsub(x = row.names(tab.1), ".1", "")
  aa <- gsub(x = aa, ".2", "")
  tab.1$cell <- aa
  # 
  # Resort factor analysis 
  tab.1$analysis <- factor(tab.1$analysis, levels = c("naive", "grouped_pr", "grouped_pd"))
  # Rename percentile columns 
  colnames(tab.1)[1] <- "lower"
  colnames(tab.1)[3] <- "upper"
  # 
  p<- ggplot(tab.1, aes(x=analysis, y=Estimate, group=cell)) +
    geom_line(aes(color=cell))+
    geom_point(aes(color=cell)) + facet_grid(cols =  vars(cell)) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=0.5), legend.position="bottom") + 
    xlab("") + geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05))
  p1<- ggplot(tab.1, aes(x=analysis, y=Estimate, group=cell)) +
    geom_bar(stat = "identity", aes(fill=cell)) + facet_grid(cols =  vars(cell)) + 
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=0.5), legend.position= "none") + 
    xlab("") + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05)) 
  p1
}
