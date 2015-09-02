library(ggplot2)
require(reshape2)
require(data.table)
require(limma)
library(doMC) 
registerDoMC(3) 

# read raw counts and sample info
GTEx.expr <- fread('GTEx_Analysis_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct')
sample.info <- fread('GTEx_Data_2014-01-17_Annotations_SampleAttributesDS.csv') 
sample.info <- sample.info[, .(SAMPID,SMTS,SMTSD)]

# get only data for heart left ventricle
col.list <- sample.info[SMTSD=="Heart - Left Ventricle",SAMPID]
expr <- GTEx.expr[, intersect(names(GTEx.expr), col.list), with=FALSE]

# normalize the data using voom
expr <- voom(expr, plot = T, normalize.method = "quantile")
expr <- as.data.table(expr$E)

# add expression, gene id and symbol
GTEx.expr.LV <- cbind(GTEx.expr[,.(Name,Description)],expr)
GTEx.expr.LV$Name <- sub("[.][0-9]+","",GTEx.expr.LV$Name)

# transform the expression data into a matrix 
expr.m <- melt(GTEx.expr.LV,variable.name = 'SUBJID')
expr.m <- dcast(expr.m,SUBJID~Name,value.var='value')
expr.m$SUBJID <- sub("[-][0-9]+[-].*","",expr.m$SUBJID)

# get subject info
subj.info <- fread("GTEx_Analysis_Annotations_Subject_DS__Pilot_2013_01_31.txt")
subj.info <- subj.info[,.(SUBJID,GENDER,AGE)]
subj.info$AGE <- sub(" years","",subj.info$AGE)
subj.info <- subj.info[!AGE==""]
subj.info <- cbind(subj.info,colsplit(subj.info$AGE,"-",c("min","max")))
subj.info$meanAGE <- round((subj.info$min+subj.info$max)/2)

# merge subject info to expression
expr.m <- merge(subj.info[, .(SUBJID,GENDER,meanAGE)],expr.m,by='SUBJID')

# This is not required right now
# but this step is important to reduce the size of the dataset
# filter out certain low expressing genes
# expr.m <- expr.m[,4:ncol(expr.m)]
# expr.m <- expr.m[,which(colMeans(expr.m > -3.287463) > 0.8)]
# expr.m <- cbind(expr.m[,c(1:3)],expr.m)

# number of genes to test for joint effects
len <- length(grep('ENSG',colnames(expr.m)))

# for e.g. we want to test for joint effects of MTSS1 & LINC00964 with all other genes in the genome
# Their id is: ENSG00000249816 and ENSG00000170873

# we make a 2 column table of Genes and Effector (it is called lnc here)
table <- data.frame(lnc = c(rep('ENSG00000249816',len), rep('ENSG00000170873',len)), 
                    Genes = rep(grep('ENSG',colnames(expr.m),value=T),2))

alldata <- expr.m
setkey(alldata,'SUBJID')

# Joint Effects Test function
DiffTest <- function(x,expr.sub)
{
  expr = paste("a=alldata[,list(SUBJID,GENDER,meanAGE,",as.character(x$Genes[1]),',',as.character(x$lnc[1]),")]",sep='')
  eval(parse(text=expr))
  setnames(a,4,'signal')
  setnames(a,5,'lnc')
  
  mod1a = with(a, lm(signal ~ meanAGE + GENDER + lnc))
  l1 = summary(mod1a)
  
  Corr = cor(x = a[,signal], y = a[,lnc], method = "spearman")
  meanGene = mean(a[,signal])
  meanLnc = mean(a[,lnc])
  Gene = as.character(expr.sub[which(expr.sub$Name %in% as.character(x$Genes[1])),2])
  Lnc = as.character(expr.sub[which(expr.sub$Name %in% as.character(x$lnc[1])),2])
  
  res = data.frame(
    Gene = Gene,
    Lnc = Lnc,
    SpearmanCorr = round(Corr,3),
    Mean_Gene = meanGene,
    Mean_Lnc = meanLnc,
    Lnc_Estimate = coef(l1)[4,1],
    Lnc_SE = coef(l1)[4,2],
    P = coef(l1)[4,4]
  )
  
  if(res$P[1]<1e-09)
  {
    jpeg(file=paste(Gene,'_',Lnc,'.jpg',sep=''),800,800)
    p = ggplot(a, aes(x = signal,y = lnc))+
      geom_point(size = 4)+
      xlab(paste("\n",Gene))+
      ylab(paste(Lnc,"\n"))+
      geom_smooth(method = lm)+
      ggtitle(paste("P=",res$P[1],"\nSpearman=",round(Corr,3)))+
      theme(axis.title.x=element_text(size=14,color="black"),
            axis.text.x=element_text(size=14,color="black"),
            axis.title.y=element_text(size=14,color="black"),
            axis.text.y=element_text(size=14,color="black"),
            plot.title=element_text(size=16),
            legend.title=element_text(size=12),
            legend.text=element_text(size=12))
    print(p)
    dev.off()
  } 
  
  return(res)
}

# Get name and description
expr.sub <- as.data.frame(GTEx.expr.LV[,.(Name,Description)])

# run function in parallel mode
results.JointEffects <- ddply(table,.(lnc,Genes),.fun = function(x) DiffTest(x,expr.sub),.parallel=TRUE)
results.JointEffects$P.Adj <- p.adjust(results.JointEffects$P, method='fdr')

