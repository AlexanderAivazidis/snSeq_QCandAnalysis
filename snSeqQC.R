### snSeqQC
require(Matrix)
library(Seurat)
require(myUtils)
require(Biomart)
require(ggplot2)
require(Hmisc)
library(SoupX)
require(biomaRt)

## Load smartSeq data:
data0 = read.delim('/home/jovyan/data/snSeq/smartSeq/combined/study5705-tic281-salmon-genecounts.txt', header = TRUE, row.names = 1)

# Load snSeq data:
data1 = readMM('data/snSeqQC/cellranger302_count_29507_5705STDY7945423_mm10-3_0_0_premrna/filtered_feature_bc_matrix/matrix.mtx.gz')
rowdata1 = read.delim('data/snSeqQC/cellranger302_count_29507_5705STDY7945423_mm10-3_0_0_premrna/filtered_feature_bc_matrix/features.tsv.gz', header = FALSE)
coldata1 = read.delim('data/snSeqQC/cellranger302_count_29507_5705STDY7945423_mm10-3_0_0_premrna/filtered_feature_bc_matrix/barcodes.tsv.gz', header = FALSE)
data1 = as.matrix(data1)
rownames(data1) = rowdata1[,2]
colnames(data1) = coldata1[,1]

# Load snSeq data:
data2 = readMM('data/snSeqQC/matrix.mtx.gz')
rowdata2 = read.delim('data/snSeqQC/features.tsv.gz', header = FALSE)
coldata2 = read.delim('data/snSeqQC/barcodes.tsv.gz', header = FALSE)
data2 = as.matrix(data2)
rownames(data2) = rowdata2[,2]
colnames(data2) = coldata2[,1]
colnames(data2) = paste(colnames(data2),'2', sep = '')

# Mitochondrial Genes:
ensembl = useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
mtGenes = getBM(attributes = 'mgi_symbol', filter = 'chromosome_name', values = "MT", mart = ensembl)
mtGenes = as.character(unlist(mtGenes))[as.character(unlist(mtGenes)) %in% rownames(data1)]
mtGenes_eg = c('ENSMUSG00000064341', 'ENSMUSG00000064345', 'ENSMUSG00000064351', 'ENSMUSG00000064354', 'ENSMUSG00000064356',
                   'ENSMUSG00000064357', 'ENSMUSG00000064358', 'ENSMUSG00000064360', 'ENSMUSG00000065947', 'ENSMUSG00000064363',
                   'ENSMUSG00000064367', 'ENSMUSG00000064368', 'ENSMUSG00000064370')
mtGenes_et =  c('ENSMUST00000082392.1', 'ENSMUST00000082396.1', 'ENSMUST00000082402.1', 'ENSMUST00000082405.1', 'ENSMUST00000082407.1',
                'ENSMUST00000082408.1', 'ENSMUST00000082409.1', 'ENSMUST00000082411.1', 'ENSMUST00000084013.1', 'ENSMUST00000082414.1',
                'ENSMUST00000082418.1', 'ENSMUST00000082419.1', 'ENSMUST00000082421.1')
mtProp1 = colSums(data1[mtGenes,])/colSums(data1)
mtProp2 = colSums(data2[mtGenes,])/colSums(data2)
mtProp3 = colSums(data0[mtGenes_et,])/colSums(data0)

mtDataFrame = data.frame(mtProp = c(mtProp1, mtProp2, mtProp3),
                         Run = c(rep('FFT4G_10x', length(mtProp1)), rep('OCT1_10x', length(mtProp2)), rep('FFT4G_SmartSeq', length(mtProp3))))

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
  annotate("text", x = 1.85, y = 0.00001, label = "n(y = 0) = 1284")
pdf(file = 'mitochondrialPercentage.pdf', width = 10, height = 10)
p0
dev.off()

# Now Load gene counts from STAR:
data0 = read.delim('/home/jovyan/data/snSeq/smartSeq/combined/study5705-tic281-star-genecounts.txt', header = TRUE, row.names = 1)
rownames(data0) = mapIdsMouse(rownames(data0), 'ENSEMBL', 'SYMBOL')

# Load Allen Data:
options(stringsAsFactors = FALSE)
dataAllen = read.delim('data/Allen/mouse_VISp_2018-06-14_exon+intron_cpm-matrix.csv', sep = ',')
rowdataA = read.delim('data/Allen/mouse_VISp_2018-06-14_genes-rows.csv', sep = ',')
coldataA = read.delim('data/Allen/mouse_VISp_2018-06-14_samples-columns.csv', sep = ',')
rownames(dataAllen) = rowdataA[,1]
commonGenes = intersect(intersect(rownames(dataAllen), rownames(data1)), rownames(data2))
celltypes = coldataA[,'class']
keep = celltypes %in% c('GABAergic', 'Endothelial', 'Glutamatergic', 'Non-Neuronal')
coldataA = coldataA[keep,]
dataAllen = dataAllen[,keep]
celltypes = celltypes[keep]
subtypes = unlist(lapply(1:dim(coldataA)[1], function(x) paste(coldataA[x,c('class', 'subclass')], sep = '_', collapse = '_')))

data1 = data1[commonGenes,]
data2 = data2[commonGenes,]
dataAllen = dataAllen[commonGenes,]
dataAllen = dataAllen[,2:dim(dataAllen)[2]]

# Quality Control Statistics:

genesDetected1 = colSums(data1 > 0)
genesDetected2 = colSums(data2 > 0)

readsDetected1 = colSums(data1)
readsDetected2 = colSums(data2)

numberOfCells1 = length(genesDetected1)
numberOfCells2 = length(genesDetected2)

qcDataFrame = data.frame(genesDetected = c(genesDetected1, genesDetected2),
                         readsDetected = c(readsDetected1, readsDetected2),
                         Run = c(rep('FFT4G', length(genesDetected1)), rep('OCT1', length(genesDetected2))))

p1 <- ggplot(qcDataFrame, aes(x=Run, y=genesDetected, fill = Run)) + 
  geom_violin() + 
geom_jitter(shape=16, size = 0.25, width = 0.45, height = 0.45) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14)) +
  annotate("text", x = 0.7, y = 5000, label = "n = 4865") +
  annotate("text", x = 1.85, y = 5000, label = "n = 5814") + ylab('Genes Detected')
pdf(file = 'NumberOfDetectedGenes.pdf', width = 10, height = 10)
p1
dev.off()

# Load number of doublets from Srublet run:

doubletRate1 = as.double(read.table('/home/jovyan/data/snSeq/FFT4G/filtered_feature_bc_matrix/overall_doublets_rate.txt'))
doubletRate2 = as.double(read.table('/home/jovyan/data/snSeq/OCT1/cellranger302_count_29507_5705STDY7945424_mm10-3_0_0_premrna/filtered_feature_bc_matrix/overall_doublets_rate.txt'))

df <- data.frame(Run=c("FFT4G", "OCT1"),
                 Doublet_Rate=c(doubletRate1, doubletRate2))|
pA <-ggplot(data=df, aes(x=Run, y=Doublet_Rate, fill = Run)) + ylab('Doublet Rate') +
  geom_bar(stat="identity") +
theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"),
      legend.text=element_text(size=14))
pdf(file = 'DoubletRate.pdf', width = 10, height = 10)
pA
dev.off()

### Remove soup: ### ?????

# Load snSeq data as soupX object:

scl_FFT4G = load10X('/home/jovyan/data/snSeq/FFT4G/')
scl_OCT1 = load10X('/home/jovyan/data/snSeq/OCT1/cellranger302_count_29507_5705STDY7945424_mm10-3_0_0_premrna/')

# For each cluster identify genes that should not be expressed in this cluster 
# but are highly expressed in other clusters, based on the Allen data set:

gaba = c('Slc6a1', 'Slc6a13', 'Slc6a11')
glutamat = c('Slc17a7', 'Slc17a6', 'Grin1')

classes = unique(coldataA[,'class'])
classes = classes[c(1,3)]
labs = list()
for (i in 1:length(classes)){
  labs[[i]] = colnames(scl_FFT4G$toc)[predictions1.1 == classes[i]]
}
toUse = matrix(unlist(lapply(1:length(classes), function(x) colnames(scl_FFT4G$toc) %in% labs[[x]])), nrow = length(classes),
               dimnames = list(classes, colnames(scl_FFT4G$toc)), byrow = TRUE)
unexpressed = list("GABAergic" = glutamat, "Endothelial" = c(glutamat, gaba), "Glutamatergic" = gaba, "Non-Neuronal" = c(glutamat, gaba))

scl_FFT4G = calculateContaminationFraction(scl_FFT4G, "Channel1", list('glut' = glutamat, 'gaba' = gaba))
scl_OCT1 = calculateContaminationFraction(scl_OCT1, "Channel1", list('glut' = glutamat, 'gaba' = gaba))

plotChannelContamination(scl_FFT4G, "Channel1")
plotChannelContamination(scl_OCT1, "Channel1")

scl_FFT4G = interpolateCellContamination(scl_FFT4G, "Channel1")
scl_OCT1 = interpolateCellContamination(scl_OCT1, "Channel1")

scl_FFT4G = adjustCounts(scl_FFT4G)
scl_OCT1 = adjustCounts(scl_OCT1)

# Remove doublet cells:

doublets_FFT4G = read.table('/home/jovyan/data/snSeq/FFT4G/filtered_feature_bc_matrix/predicted_doublets.txt')
doublets_OCT1 = read.table('/home/jovyan/data/snSeq/OCT1/cellranger302_count_29507_5705STDY7945424_mm10-3_0_0_premrna/filtered_feature_bc_matrix/predicted_doublets.txt')

data1 = data1[, doublets_FFT4G == 0]
data2 = data2[, doublets_OCT1 == 0]

# Classify cell types:

visualCortexData = cbind(dataAllen,data1,data2)
metaData = cbind(c(rep('Allen', dim(dataAllen)[2]), rep('snSeq1', dim(data1)[2]), rep('snSeq2', dim(data2)[2])),
                 c(coldataA[,'class'], rep('Unknown', dim(data1)[2]), rep('Unknown', dim(data2)[2])),
                 c(subtypes, rep('Unknown', dim(data1)[2]), rep('Unknown', dim(data2)[2])))
colnames(metaData) = c('tech', 'celltype', 'subtype')
rownames(metaData) = c(colnames(dataAllen), colnames(data1), colnames(data2))
metaData = as.data.frame(metaData)

visualCortex <- CreateSeuratObject(visualCortexData, meta.data = metaData)
visualCortex.list <- SplitObject(visualCortex, split.by = 'tech')

for (i in 1:length(visualCortex.list)) {
  visualCortex.list[[i]] <- NormalizeData(visualCortex.list[[i]], verbose = FALSE)
  visualCortex.list[[i]] <- FindVariableFeatures(visualCortex.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                                 verbose = FALSE)
}

visualCortex.query <- visualCortex.list[["snSeq1"]]
visualCortex.anchors <- FindTransferAnchors(reference = visualCortex.list[['Allen']], query = visualCortex.query, 
                                            dims = 1:30)
predictions <- TransferData(anchorset = visualCortex.anchors, refdata = visualCortex.list[['Allen']]$celltype, 
                            dims = 1:30)
predictions1.1 = predictions[,1]
predictions_Subtype <- TransferData(anchorset = visualCortex.anchors, refdata = visualCortex.list[['Allen']]$subtype, 
                                    dims = 1:30)
predictions1.2 = predictions_Subtype[,1]

scores1.1 = colSums(predictions[,2:5])
scores1.1 = scores1.1/sum(scores1.1)
names(scores1.1) = c('GABAergic', 'Endothelial','Glutamatergic', 'NonNeuronal')
scores1.1 = scores1.1[order(scores1.1)]
saveRDS(scores1.1, file = 'scores1.1.rds')

scores1.2 = colSums(predictions_Subtype[,2:24])
scores1.2 = scores1.2/sum(scores1.2)
names(scores1.2) = unlist(lapply(1:length(scores1.2), function(x) paste(strsplit(names(scores1.2)[x], split = '\\.')[[1]][-c(1,2)], sep = '', collapse = '')))
scores1.2 = scores1.2[order(names(scores1.2))]
saveRDS(scores1.2, file = 'scores1.2.rds')

visualCortex.query <- visualCortex.list[["snSeq2"]]
visualCortex.anchors <- FindTransferAnchors(reference = visualCortex.list[['Allen']], query = visualCortex.query, 
                                            dims = 1:30)
predictions <- TransferData(anchorset = visualCortex.anchors, refdata = visualCortex.list[['Allen']]$celltype, 
                            dims = 1:30)
predictions2.1 = predictions[,1]
predictions_Subtype <- TransferData(anchorset = visualCortex.anchors, refdata = visualCortex.list[['Allen']]$subtype, 
                                    dims = 1:30)
predictions2.2 = predictions_Subtype[,1]

scores2.1 = colSums(predictions[,2:5])
scores2.1 = scores2.1/sum(scores2.1)
names(scores2.1) = c('GABAergic', 'Endothelial','Glutamatergic', 'NonNeuronal')
scores2.1 = scores2.1[order(scores2.1)]
saveRDS(scores2.1, file = 'scores2.1.rds')

scores2.2 = colSums(predictions_Subtype[,2:24])
scores2.2 = scores2.2/sum(scores2.2)
names(scores2.2) = unlist(lapply(1:length(scores2.2), function(x) paste(strsplit(names(scores2.2)[x], split = '\\.')[[1]][-c(1,2)], sep = '', collapse = '')))
scores2.2 = scores2.2[order(names(scores2.2))]
saveRDS(scores2.2, file = 'scores2.2.rds')

scoresAllen = table(coldataA['class'])
scoresAllen = scoresAllen/sum(scoresAllen)
names(scoresAllen) = c('Endothelial', 'GABAergic', 'Glutamatergic', 'NonNeuronal')
scoresAllen = scoresAllen[names(scores1.1)]

subtypes = unlist(lapply(1:length(scores1.2), function(x) paste(strsplit(names(scores1.2)[x], split = '\\.')[[1]][-c(1,2)], sep = '', collapse = '')))
scoresAllen.2 = table(subtypes)
scoresAllen.2 = scoresAllen.2/sum(scoresAllen.2)
scoresAllen.2 = scoresAllen.2[sort(names(scoresAllen.2))]
missing1 = names(scoresAllen.2)[which(!names(scoresAllen.2) %in% names(scores2.2))] 
missing2 = names(scores2.2)[!names(scores2.2) %in% names(scoresAllen.2)] # I checked these are in the same order
for (i in 1:length(missing1)){
  names(scoresAllen.2)[names(scoresAllen.2) == missing1[i]] = missing2[i]
}
scoresAllen.2 = scoresAllen.2[names(scores1.2)]

allScores.1 = data.frame(Celltype = names(scores1.1), Proportion = c(scoresAllen, scores1.1, scores2.1), Run = c(rep('Allen', length(scoresAllen)),
                                                                                                                 rep('FFT4G', length(scores1.1)), rep('OCT1', length(scores2.1))))

allScores.2 = data.frame(Celltype = names(scores1.2), Proportion = c(scoresAllen.2, scores1.2, scores2.2), Run = c(rep('Allen', length(scoresAllen.2)),
                                                                                                                   rep('FFT4G', length(scores1.2)), rep('OCT1', length(scores2.2))))

library(tidyr)
library(ggplot2)

p2 = ggplot(data = allScores.1, 
       aes(x = Celltype, y = Proportion, fill = Run)) + 
  geom_bar(stat = 'identity', position = 'dodge')

p3 = ggplot(data = allScores.2, 
       aes(x = Celltype, y = Proportion, fill = Run)) + 
  geom_bar(stat = 'identity', position = 'dodge') + theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf(file = 'Classification_CellClasses.pdf', width = 10, height = 5)
p2
dev.off()

pdf(file = 'Classification_CellSubClasses.pdf', width = 10, height = 5)
p3
dev.off()

visualCortex.list[[2]]$celltype = predictions1.1
visualCortex.list[[3]]$celltype = predictions2.1

visualCortex.list[[2]]$subtype = predictions1.2
visualCortex.list[[3]]$subtype = predictions2.2

## Now do clustering with all three datasets
reference.list <- visualCortex.list[c("Allen", "snSeq1", "snSeq2")]
visualCortex.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
visualCortex.integrated <- IntegrateData(anchorset = visualCortex.anchors, dims = 1:30)

library(ggplot2)
library(cowplot)
# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(visualCortex.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
visualCortex.integrated <- ScaleData(visualCortex.integrated, verbose = FALSE)
visualCortex.integrated <- RunPCA(visualCortex.integrated, npcs = 30, verbose = FALSE)
visualCortex.integrated <- RunUMAP(visualCortex.integrated, reduction = "pca", dims = 1:30)
visualCortex.integrated$tech[visualCortex.integrated$tech == 'snSeq1'] = 'FFT4G'
visualCortex.integrated$tech[visualCortex.integrated$tech == 'snSeq2'] = 'OCT1'
# VisualCortex.snSeq = subset(visualCortex.integrated, cells = c(colnames(data1), colnames(data2)))
p4 <- DimPlot(visualCortex.integrated, reduction = "umap", group.by = "tech")
p5 <- DimPlot(visualCortex.integrated, reduction = "umap", group.by = "celltype")
p6 <- DimPlot(visualCortex.integrated, reduction = "umap", group.by = "subtype")
p15 = FeaturePlot(visualCortex.integrated, feature = "Rorb", reduction = "umap", cols = c('yellow', 'red'))
p19 = FeaturePlot(visualCortex.integrated, feature = "nFeature_RNA", reduction = "umap", cols = c('yellow', 'red')) + labs(title =   'Number of Genes Detected')

pplot_grid(p4, p6, p5, align = 'v')

pdf(file = 'Clustering_Run.pdf', width = 7, height = 7)
p4
dev.off()

pdf(file = 'Clustering_CellClass.pdf', width = 7, height = 7)
p5
dev.off()

pdf(file = 'Clustering_CellSubClass.pdf', width = 12, height = 7)
p6
dev.off()

pdf(file = 'Clustering_Rorb.pdf', width = 7, height = 7)
p15
dev.off()

pdf(file = 'Clustering_NumberOfDetectedGenes.pdf', width = 7, height = 7)
p19
dev.off()

### Make a second set of plots with Allen data excluded:
p7 <- DimPlot(visualCortex.integrated[,visualCortex.integrated$tech %in% c('FFT4G', 'OCT1')], reduction = "umap", group.by = "tech", pt.size = 0.05)
p8 <- DimPlot(visualCortex.integrated[,visualCortex.integrated$tech %in% c('FFT4G', 'OCT1')], reduction = "umap", group.by = "celltype", pt.size = 0.05)
p9 <- DimPlot(visualCortex.integrated[,visualCortex.integrated$tech %in% c('FFT4G', 'OCT1')], reduction = "umap", group.by = "subtype", pt.size = 0.05)
p13 = FeaturePlot(visualCortex.integrated, feature = "Rorb", reduction = "umap")
### Make a third set of plots with Allen data excluded from common space:

## Now do clustering with all three datasets
reference.list <- visualCortex.list[c("snSeq1", "snSeq2")]
visualCortex.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
visualCortex.integrated <- IntegrateData(anchorset = visualCortex.anchors, dims = 1:30)
# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(visualCortex.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
visualCortex.integrated <- ScaleData(visualCortex.integrated, verbose = FALSE)
visualCortex.integrated <- RunPCA(visualCortex.integrated, npcs = 30, verbose = FALSE)
visualCortex.integrated <- RunUMAP(visualCortex.integrated, reduction = "pca", dims = 1:30)
visualCortex.integrated$tech[visualCortex.integrated$tech == 'snSeq1'] = 'FFT4G'
visualCortex.integrated$tech[visualCortex.integrated$tech == 'snSeq2'] = 'OCT1'
# VisualCortex.snSeq = subset(visualCortex.integrated, cells = c(colnames(data1), colnames(data2)))
p10 <- DimPlot(visualCortex.integrated, reduction = "umap", group.by = "tech")
p11 <- DimPlot(visualCortex.integrated, reduction = "umap", group.by = "celltype")
p12 <- DimPlot(visualCortex.integrated, reduction = "umap", group.by = "subtype")
p13 = FeaturePlot(visualCortex.integrated, feature = "Rorb", reduction = "umap")
p14 = FeaturePlot(visualCortex.integrated, feature = "Rorb", reduction = "umap", cols = c('yellow', 'red'))
plot_grid(p4, p6, p5, align = 'v')

pdf(file = 'Clustering_Run_snSeqOnly.pdf', width = 7, height = 7)
p10
dev.off()

pdf(file = 'Clustering_CellClass_snSeqOnly.pdf', width = 7, height = 7)
p11
dev.off()

pdf(file = 'Clustering_CellSubClass_snSeqOnly.pdf', width = 10, height = 7)
p12
dev.off()

pdf(file = 'Clustering_RorbExpression_snSeqOnly.pdf', width = 7, height = 7)
p13
dev.off()

pdf(file = 'Clustering_RorbExpression_snSeqOnly_yellowToRed.pdf', width = 7, height = 7)
p14
dev.off()


