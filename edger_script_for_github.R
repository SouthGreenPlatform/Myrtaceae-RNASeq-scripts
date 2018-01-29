# EdgeR based pipeline used for RNAseq data analysis of the response 
# of myrtaceae species (close species to Eucalyptus grandis) to a fungi attack

# R packages to load 
library(QuasR)
library(GenomicFeatures)
library(corrplot)
library(ggbio)
library(edgeR)
library(lattice)
library(gplots)
library(GenomicAlignments)
library(reshape2)

# Set the wording directory
setwd("...")

# Input a tabulated file of the counts, having the structure below 
# gene gene_width  sample1  sample2  sample3  sample4 sample5 ...

gene         width   V3   V4   V1   V2   V6   V7
1 LOC104414021  2543 1827 2158 1404 1845 2077 2126
2 LOC104414022  2924  144  147  116  218  105  165
3 LOC104414023  2234  213  148  136  137  212  220
...

table_of_counts <- read.delim("counts_trigl_bds_eucgr.lst")
names(table_of_counts)
[1] "gene"  "width" "V3"    "V4"    "V1"    "V2"    "V6"    "V7"

# RPKM calculation
rpkmDFgene <- t(t(table_of_counts[,3:8]/table_of_counts[,2] * 1000)/colSums(table_of_counts[,3:8]) *1e6)
rownames(rpkmDFgene) <- table_of_counts[,1]

V3        V4        V1        V2        V6        V7
LOC104414021 45.677043 57.025195 51.315046 38.457457 41.323088 48.994153
LOC104414022  3.131057  3.378327  3.687267  3.951934  1.816831  3.306998
LOC104414023  6.061809  4.451847  5.658219  3.250633  4.801261  5.771210
...

write.table(rpkmDFgene, "trigl_bds_eucgr_rpkmDF.lst", quote=FALSE, sep="\t", col.names = NA)

countDF <- table_of_counts[,3:8]
rownames(countDF) <- table_of_counts[, 1]

# Calculation of the genewise standard deviation for the counts
countDF=cbind(countDF,apply(countDF[,1:6],1,sd))

# Calculation of the genewise mean for the counts
countDF=cbind(countDF,apply(countDF[,1:6],1,mean)) 

# Plot mean-sd relationship
plot(countDF[,8],countDF[,7])

# Plot the correlations between samples
d <- cor(rpkmDFgene, method="spearman")
corrplot(d)

# Plot a dendrogram of the samples to see their clustering
hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)

# Experimental design  
# 1 factor : reaction to fungi attack
# 2 levels : tolerant(resistant), sensible
# 6 samples : 2 tolerant - 6 sensible

# Indicate the sample status. The column order is the same as in table_of_counts
group <- c("tolerant","tolerant","sensible","sensible","sensible","sensible")

y <- DGEList(counts=countDF[,1:6],group=group)

# Keep only genes with cmp()>1 and expressed in at least in 2 samples 
keep=rowSums(cpm(y)>1) >=2 
y=y[keep,]


# Multidimensional scaling plot of distances between digital gene expression profiles
# The distances on the plot approximate the typical log2 fold changes between the samples.
mds <- plotMDS(y, main="MultiDimensional Scaling plot", col=as.numeric(y$samples$group), xlab="Axe1 : logFC", ylab="Axe2 : logFC")
legend("topright",levels(y$samples$group),col=1:4,pch=20,cex=0.6)

# Plot the log-fold change (i.e. the log of the ratio of expression levels for each gene between two experimential groups) 
# against the log-concentration (i.e. the overall average expression level for each gene across the two groups). 
plotSmear(y,xlab="Average logCPM", ylab="logFC")

# Normalisation of the libraries
y<-calcNormFactors(y)

# Estimation of the common dispersion
y <- estimateCommonDisp(y)

# Estimation of the tagwise dispersion
# y <- estimateTagwiseDisp(y)	## To calculate a specific dispersion for each tag instead of using the same dispersion for all the tags
				## the common dispersion must be always calculate before the tagwise dispersion

# Statistical tests based on negative binomial distribution
et <- exactTest(y, pair=c("tolerant", "sensible"), dispersion="common")

# Calculation of the number of under or over expressed genes 
# treshold p.value=0.10 and absolute value of logFC=1. 
# Due to multiple testing p-values are adjusted with FDR method
summary(de <- decideTestsDGE(et,adjust.method="fdr",p.value=0.10, lfc=1))

[,1] 
-1   125
0  42690
1    617

number of DE genes : 125 + 617 = 742

# Extracts the genes ranked by p-value 
edge <- as.data.frame(topTags(et, n=43432))

## Select genes with logFC >= 1 or logFC <= -1
edge2fold <- edge[edge$logFC >= 1 | edge$logFC <= -1,]

## then select genes with FDR < 0.10
edge2foldpadj <- edge2fold[edge2fold$FDR <= 0.10, ]
nrow(edge2foldpadj)
[1] 742

# Write DE genes in a file
write.table(edge2foldpadj, "trista_bds_eucgr_DE_genes_datapaper.lst", row.names=FALSE, quote=FALSE, sep="\t")

# Plot a heatmap 
rpkm_scaled <- rpkmDFgene[rownames(edge2foldpadj)[1:742],c(1,2,3,4,5,6)] 
rpkm_scaled <- rpkm_scaled[order(rpkm_scaled[,1]),]

colnames(rpkm_scaled) <- c("tolerant","tolerant","sensible","sensible","sensible","sensible")
rpkm_scaled <- t(scale(t(as.matrix(rpkm_scaled))))

levelplot(t(rpkm_scaled[1:30,]), height=0.2, col.regions=colorpanel(40, "darkblue", "yellow", "white"), main="Expression Values (DEG Filter: FDR 1%, FC > 2)", colorkey=list(space="top"), xlab="", ylab="Gene ID")

# Plot the expression level by conditions for the top most DE genes
# the first 30 most DE genes will be displayed
y_long= melt(cpm(y[,c(1,2,3,4,5,6)])) 
colnames(y_long) <- c("gene","condition","rpkm")
table(y_long$condition)

y_top_de = y_long[y_long$gene %in% rownames(edge2foldpadj)[1:30],]

ggplot(y_top_de,aes(x=condition,y=rpkm,colour=gene,group=gene))+geom_point()+geom_smooth(method="lm",se=F)+scale_y_log10()

# Load the annotation file of the Eucalyptus reference genome
txdb <- makeTxDbFromGFF(file="/bank/genfam/genome_data_v1/REF_EUCGR/EUCGR-NCBI100-sequence_feature-genfam.gff3",format="gff3",dataSource="NBBI100",organism="Eucalyptus grandis")
saveDb(txdb, file="/work/NGSwaiting4newNAS/myrtaceae/R_DATA/syzygium/edger/ezucgr.sqlite")

# List the columns of the annotation file
columns(txdb)
[1] "CDSCHROM"   "CDSEND"     "CDSID"      "CDSNAME"    "CDSSTART"   "CDSSTRAND" 
[7] "EXONCHROM"  "EXONEND"    "EXONID"     "EXONNAME"   "EXONRANK"   "EXONSTART" 
[13] "EXONSTRAND" "GENEID"     "TXCHROM"    "TXEND"      "TXID"       "TXNAME"    
[19] "TXSTART"    "TXSTRAND"   "TXTYPE"    

# Select the DE genes
most_de_id=rownames(edge2foldpadj[1:742,])

# List the DE genes annotations
most_de_loc=select(txdb,most_de_id,columns=c("GENEID","TXNAME","TXCHROM","TXSTART","TXEND"),keytype="GENEID")













