import numpy as np
import sys
import os
from sklearn import ensemble

if sys.platform == "darwin":
  base = "/Volumes/cycle/rrna/"
else:
  base = "/homes/gws/" + os.getlogin() + "/rrna/" # make portable for katies and svens

source_home = base + "src/"
data_home = base + "data/"

genes = np.load(data_home + "featureExtractors/Homo_sapiens.GRCh38.84_ordered_gene_names.npy")
rrna_genes = np.load(data_home + "/annotations/Homo_sapiens.GRCh38.84_rrna_gene_names.npy")
mrp = [name for name in genes if name.startswith("MRP")]
rp = [name for name in genes if name.startswith("RP")]
rrna_genes = np.hstack([rrna_genes, mrp, rp])


finished_genes = [os.path.splitext(i)[0] for i in os.listdir(data_home + "/rnafold")]
finished_genes_rrna = np.zeros(len(finished_genes))

row_names = np.load(data_home + "rnafold_results/row_names.npy")
print len(row_names)
mfe = np.load(data_home + "rnafold_results/mfe.npy")
pair_probs = np.load(data_home + "rnafold_results/pair_probs.npy")

rrna_labels = np.zeros(len(row_names))
for i,gene in enumerate(row_names):
  if len(gene) > 1:
    if gene[0] in rrna_genes:
        rrna_labels[i] = 1
    else:
        rrna_labels[i] = 0
  # else: 
    # print gene
        
np.save(data_home + "/rnafold_results/rrna_labels.npy",finished_genes_rrna)


# Predict pair probabilities
def train_and_predict_gb():
  gb = ensemble.GradientBoostingRegressor(n_estimators=200, max_depth=3) 
  gb.fit(pair_probs, rrna_labels)
  gb_pred = gb.predict(pair_probs)
  rms = rmse(gb_pred, rrna_labels)

