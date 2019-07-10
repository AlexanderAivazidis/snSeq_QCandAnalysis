# snSeq_QCandAnalysis
Standard pipeline for downstream analysis of single nucleus sequencing data

To run on Jupyter hub first copy paste the following commands into the terminal to get data from our /nfs storage on the farm:

cd
mkdir data/snQCandAnalysis/
mkdir data/snQCandAnalysis/OCT1_10x/
mkdir data/snQCandAnalysis/FFT4G_10x/
mkdir data/snQCandAnalysis/FFT4G_NEB/

Then replace aa16 in the following commands with your own sangerID and type in your farm3 password to access the sequencing data:

rsync -avzhr aa16@farm3-login:/nfs/team283/sequencing/snSeq/cellranger302_count_29507_5705STDY7945423_mm10-3_0_0_premrna/filtered_feature_bc_matrix data/snQCandAnalysis/FFT4G_10x/
gunzip data/snQCandAnalysis/FFT4G_10x/filtered_feature_bc_matrix/features.tsv.gz

rsync -avzhr aa16@farm3-login:/nfs/team283/sequencing/snSeq/cellranger302_count_29507_5705STDY7945424_mm10-3_0_0_premrna/filtered_feature_bc_matrix data/snQCandAnalysis/OCT1_10x/
gunzip data/snQCandAnalysis/OCT1_10x/filtered_feature_bc_matrix/features.tsv.gz

rsync -avzhr aa16@farm3-login:/nfs/team283/sequencing/snSeq/tic-281/combined_premrna/premrna-study5705-tic281-star-genecounts.txt data/snQCandAnalysis/FFT4G_NEB/

Finally clone the repository:

chmod 600 ~/.ssh/id_rsa
git clone git@github.com:AlexanderAivazidis/snSeq_QCandAnalysis.git
cd snSeq_QCandAnalysis/
