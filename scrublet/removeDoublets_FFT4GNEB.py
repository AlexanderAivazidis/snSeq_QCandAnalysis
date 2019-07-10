import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

tag = 'FFT4G_NEB'
output_dir = '/home/jovyan/snSeq_QCandAnalysis/scrublet'
input_dir = '/home/jovyan/data/snQCandAnalysis/FFT4G_NEB/'
input_file = '/home/jovyan/data/snQCandAnalysis/FFT4G_NEB/premrna-study5705-tic281-star-genecounts.txt'
counts_matrix = pd.read_csv(input_file, sep='\t')
counts_matrix = counts_matrix.iloc[0:-4,:]
genes = counts_matrix.iloc[:,0]
counts_matrix = counts_matrix.iloc[:,1:]
counts_matrix = np.asarray(counts_matrix)
counts_matrix = counts_matrix.transpose()

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06, sim_doublet_ratio = 10)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

predicted_doublets = predicted_doublets*1
predicted_doublets = predicted_doublets.astype(int)
detected_doublets_rate = round(scrub.detected_doublet_rate_, 4)
overall_doublets_rate = round(scrub.overall_doublet_rate_, 4)

np.savetxt(output_dir + '/' + tag + '_' + 'doublets_scores.txt', doublet_scores)   
np.savetxt(output_dir +  '/' + tag + '_' + 'predicted_doublets.txt', predicted_doublets)                              
with open(output_dir +  '/' + tag + '_' + 'detected_doublets_rate.txt', 'w') as f:
  f.write('%f' % detected_doublets_rate)  

with open(output_dir + '/' + tag + '_' +  'overall_doublets_rate.txt', 'w') as f:
  f.write('%f' % overall_doublets_rate)

f = scrub.plot_histogram()
f[0].savefig(output_dir + '/' + tag + '_' + "doubletScore_histogram.pdf", bbox_inches='tight')
