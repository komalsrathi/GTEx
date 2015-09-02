# prepare GTEx annotation for analysis but don't dcast (the last step) it yet

library(ggplot2)
library(reshape2)
library(data.table)

# To get scatter plots between two genes across all tissues
# get data for genes you want to make plots for
dat <- as.data.frame(dat[Description %like% "MTSS1$|LINC00964"])

# list of tissues for which the data is available
list <- as.character(unique(dat$SMTS))

# using dcast, transform the data
dat.scatter <- dcast(dat, variable + SMTS + SMTSD ~ Description, value.var="value")

# this example shows how you can make scatter plots of each tissue type and save individual files of the plots
# make scatter plots for each tissue type
for(i in 1:length(list))
{
  dat.tissue <- dat.scatter[grep(paste('^',list[i],'$',sep=""),dat.scatter$SMTS),]
  cor1 <- round(cor(x=dat.tissue$LINC00964,y=dat.tissue$MTSS1,method="pearson"),3)
  cor2 <- round(cor(x=dat.tissue$LINC00964,y=dat.tissue$MTSS1,method="spearman"),3)
  jpeg(file=paste('Tissue_specific_correlation_plots_MTSS1_LINC00964/',as.character(dat.tissue[1,2]),'_',"MTSS1",'_',"LINC00964",'.jpg',sep=''),
       800,800)
  p = ggplot(dat.tissue,aes(LINC00964,MTSS1)) + 
    theme(axis.text.x=element_text(size=12,color="black"),
          axis.text.y=element_text(size=12,color="black"),
          strip.text.x=element_text(size=18),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14)) + 
    geom_point(cex=2.5,pch=20) + facet_grid(.~SMTS) + 
    ggtitle(paste("Pearson",cor1,"\nSpearman:",cor2,"\n",sep=" "))  
  print(p)
  dev.off()
}

# this example shows how you can make boxplots of each type and save in one image
# make boxplots for each tissue type
# melt the transformed data
dat.boxplot <- melt(dat.scatter,variable.name='Description')
ggplot(dat.boxplot,aes(Description,log2(value),fill=factor(Description))) +
  theme(axis.text.x=element_text(size=12,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        strip.text.x=element_text(size=18),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14)) + 
  geom_boxplot(show_guide=F) + ggtitle("Boxplot: MTSS1 & LINC00964\n") + facet_wrap(~SMTS,ncol=5)


