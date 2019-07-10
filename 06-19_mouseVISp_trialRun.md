Single Nucleus Sequencing Trial Quality Control and Analysis
================

To reproduce this analysis on CellGenes JupyterHub copy paste the
following commands into the terminal to get data from our /nfs storage
on the farm:

cd mkdir data/snQCandAnalysis/  
mkdir data/snQCandAnalysis/OCT1\_10x/  
mkdir data/snQCandAnalysis/FFT4G\_10x/  
mkdir data/snQCandAnalysis/FFT4G\_NEB/

Then replace aa16 in the following commands with your own sangerID and
type in your farm3 password, when asked to do so, to access the
sequencing data:  
rsync -avzhr
<aa16@farm3-login>:/nfs/team283/sequencing/snSeq/cellranger302\_count\_29507\_5705STDY7945423\_mm10-3\_0\_0\_premrna/filtered\_feature\_bc\_matrix
data/snQCandAnalysis/FFT4G\_10x/  
gunzip
data/snQCandAnalysis/FFT4G\_10x/filtered\_feature\_bc\_matrix/features.tsv.gz

rsync -avzhr
<aa16@farm3-login>:/nfs/team283/sequencing/snSeq/cellranger302\_count\_29507\_5705STDY7945424\_mm10-3\_0\_0\_premrna/filtered\_feature\_bc\_matrix
data/snQCandAnalysis/OCT1\_10x/  
gunzip
data/snQCandAnalysis/OCT1\_10x/filtered\_feature\_bc\_matrix/features.tsv.gz

rsync -avzhr
<aa16@farm3-login>:/nfs/team283/sequencing/snSeq/tic-281/combined\_premrna/premrna-study5705-tic281-star-genecounts.txt
data/snQCandAnalysis/FFT4G\_NEB/

Now we load we load the three datasets into R, as well as the required R
packages:

``` r
require(Matrix)
library(Seurat)
require(myUtils)
require(ggplot2)
require(Hmisc)
library(SoupX)
require(biomaRt)

# Neb data:
data_NEB = read.delim('/home/jovyan/data/snQCandAnalysis/FFT4G_NEB/premrna-study5705-tic281-star-genecounts.txt', header = TRUE, row.names = 1)

# 10x data:
data_10xF = readMM('/home/jovyan/data/snQCandAnalysis/FFT4G_10x/filtered_feature_bc_matrix/matrix.mtx.gz')
rowdata_10xF = read.delim('/home/jovyan/data/snQCandAnalysis/FFT4G_10x/filtered_feature_bc_matrix/features.tsv', header = FALSE)
coldata_10xF = read.delim('/home/jovyan/data/snQCandAnalysis/FFT4G_10x/filtered_feature_bc_matrix/barcodes.tsv.gz', header = FALSE)
data_10xF = as.matrix(data_10xF)
rownames(data_10xF) = rowdata_10xF[,2]
colnames(data_10xF) = coldata_10xF[,1]

data_10xO = readMM('/home/jovyan/data/snQCandAnalysis/OCT1_10x/filtered_feature_bc_matrix/matrix.mtx.gz')
rowdata_10xO = read.delim('/home/jovyan/data/snQCandAnalysis/OCT1_10x/filtered_feature_bc_matrix/features.tsv', header = FALSE)
coldata_10xO = read.delim('/home/jovyan/data/snQCandAnalysis/OCT1_10x/filtered_feature_bc_matrix/barcodes.tsv.gz', header = FALSE)
data_10xO = as.matrix(data_10xO)
rownames(data_10xO) = rowdata_10xO[,2]
colnames(data_10xO) = coldata_10xO[,1]
colnames(data_10xO) = paste(colnames(data_10xO),'2', sep = '')
```

The first QC step we make is to look at the fraction of counts/reads
from mitochondrial genes and potentially remove any outliers:

``` r
# NEB
# Mitochondrial Genes:
# ensembl = useMart("ensembl")
# ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
# mtGenes = getBM(attributes = 'mgi_symbol', filter = 'chromosome_name', values = "MT", mart = ensembl)
# mtGenes = as.character(unlist(mtGenes))[as.character(unlist(mtGenes)) %in% rownames(data_10xF)]
mtGenes_symbol = c("mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2", "mt-Atp8", "mt-Atp6", "mt-Co3",
                   "mt-Nd3", "mt-Nd4l", "mt-Nd4", "mt-Nd5", "mt-Nd6", "mt-Cytb")
mtGenes_eg = c('ENSMUSG00000064341', 'ENSMUSG00000064345', 'ENSMUSG00000064351', 'ENSMUSG00000064354', 'ENSMUSG00000064356',
                   'ENSMUSG00000064357', 'ENSMUSG00000064358', 'ENSMUSG00000064360', 'ENSMUSG00000065947', 'ENSMUSG00000064363',
                   'ENSMUSG00000064367', 'ENSMUSG00000064368', 'ENSMUSG00000064370')
# mtGenes_et =  c('ENSMUST00000082392.1', 'ENSMUST00000082396.1', 'ENSMUST00000082402.1', 'ENSMUST00000082405.1', 'ENSMUST00000082407.1',
#                 'ENSMUST00000082408.1', 'ENSMUST00000082409.1', 'ENSMUST00000082411.1', 'ENSMUST00000084013.1', 'ENSMUST00000082414.1',
#                 'ENSMUST00000082418.1', 'ENSMUST00000082419.1', 'ENSMUST00000082421.1')
mtProp0 = colSums(data_NEB[mtGenes_eg,])/colSums(data_NEB)
mtProp1 = colSums(data_10xF[mtGenes_symbol,])/colSums(data_10xF)
mtProp2 = colSums(data_10xO[mtGenes_symbol,])/colSums(data_10xO)

mtDataFrame = data.frame(mtProp = c(mtProp0, mtProp1, mtProp2),
                         Run = c(rep('FFT4G_SmartSeq', length(mtProp0)), rep('FFT4G_10x', length(mtProp1)), rep('OCT1_10x', length(mtProp2))))

p0 <- ggplot(mtDataFrame, aes(x=Run, y=mtProp, fill = Run)) + 
  scale_y_log10() +
  geom_violin() + 
  geom_jitter(shape=16, size = 0.45, width = 0.45, height = 0.05) +
  stat_summary(fun.y=mean, geom="point", size=2, color="red") +
  ylab('Mitochondrial RNA Counts Proportion') +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14)) +
  annotate("text", x = 0.7, y = 0.00001, label = "n(y = 0) = 3662") +
  annotate("text", x = 1.85, y = 0.00001, label = "n(y = 0) = 160") +
  annotate("text", x = 2.85, y = 0.00001, label = "n(y = 0) = 1284")
p0
```

![](06-19_mouseVISp_trialRun_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
pdf(file = 'figures/mitochondrialPercentage.pdf', width = 10, height = 10)
p0
dev.off()
```

    ## png 
    ##   2
