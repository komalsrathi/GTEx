library(reshape2)
library(data.table)

# download the raw files from GTEx
gtex <- fread('GTEx_Analysis_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct')
gtex.melt <- melt(gtex,id.vars=c("Name","Description"),variable.name="SAMPID")
gtex.ann <- fread('GTEx_Data_2014-01-17_Annotations_SampleAttributesDS.txt')
dat <- merge(gtex.ann,gtex.melt,by="SAMPID")
dat <- dcast(dat, SAMPID + SMTS + SMTSD ~ Description, value.var="value")
