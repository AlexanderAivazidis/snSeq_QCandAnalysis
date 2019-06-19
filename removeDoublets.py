import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

input_dir = '/home/jovyan/data/snSeq/OCT1/cellranger302_count_29507_5705STDY7945424_mm10-3_0_0_premrna/filtered_feature_bc_matrix'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

predicted_doublets = predicted_doublets*1
predicted_doublets = predicted_doublets.astype(int)
detected_doublets_rate = round(scrub.detected_doublet_rate_,4)
overall_doublets_rate = round(scrub.overall_doublet_rate_,4)

np.savetxt(input_dir + '/predicted_doublets.txt', predicted_doublets)                              
with open(input_dir + '/detected_doublets_rate.txt', 'w') as f:
  f.write('%f' % detected_doublets_rate)  

with open(input_dir + '/overall_doublets_rate.txt', 'w') as f:
  f.write('%f' % overall_doublets_rate)



scrub.plot_histogram();
print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True);
    
