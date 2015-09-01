library(reshape2)
library(data.table)

# download the raw files from GTEx
gtex <- fread('GTEx_Analysis_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct') # rpkm values or
gtex <- fread('GTEx_Analysis_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct') # raw counts

# melt counts/expression values
gtex.melt <- melt(gtex,id.vars=c("Name","Description"),variable.name="SAMPID")

# get annotation
gtex.ann <- fread('GTEx_Data_2014-01-17_Annotations_SampleAttributesDS.txt')
gtex.ann <- gtex.ann[, .(SAMPID,SMTS,SMTSD)]

# merge counts/rpkm with annotation
dat <- merge(gtex.ann,gtex.melt,by="SAMPID")
dat$Name <- sub("[.][0-9]*","",dat$Name) # remove .* from ENSEMBL ID
dt.mat <- dcast(dat, SAMPID + SMTS + SMTSD ~ Name, value.var="value")

# or you can extract any tissue type and do downstream analysis
dt <-  subset(dat, SMTSD == "Whole Blood")
dt.mat <- dcast(dt, SAMPID + SMTS + SMTSD ~ Name, value.var="value") # Name is Ensembl ID
