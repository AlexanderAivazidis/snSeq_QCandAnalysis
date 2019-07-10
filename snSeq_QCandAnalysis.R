### Functions for loading data:

load10X = function(directory){
  # Load snSeq data:
  data1 = readMM('/home/jovyan/data/snSeqQC/cellranger302_count_29507_5705STDY7945423_mm10-3_0_0_premrna/filtered_feature_bc_matrix/matrix.mtx.gz')
  rowdata1 = read.delim('/home/jovyan/data/snSeqQC/cellranger302_count_29507_5705STDY7945423_mm10-3_0_0_premrna/filtered_feature_bc_matrix/features.tsv.gz', header = FALSE)
  coldata1 = read.delim('/home/jovyan/data/snSeqQC/cellranger302_count_29507_5705STDY7945423_mm10-3_0_0_premrna/filtered_feature_bc_matrix/barcodes.tsv.gz', header = FALSE)
  data1 = as.matrix(data1)
  rownames(data1) = rowdata1[,2]
  colnames(data1) = coldata1[,1]
}