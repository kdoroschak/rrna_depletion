# coding: utf-8
get_ipython().magic(u'cd rrna/src')
import kfold
import numpy as np
rrna_pair_probs = np.load("/projects/bio/rrna/data/rnafold_results/rrna_by1_pair_probs.npy") # these will change for mouse
not_rrna_pair_probs = np.load("/projects/bio/rrna/data/rnafold_results/not_rrna_mean_pair_probs.npy")
not_rrna_pair_probs.shape
rrna_pair_probs.shape
kfold.kfold(rrna_pair_probs, not_rrna_pair_probs, save_folder="/projects/bio/rrna/data/rnafold_results/", n_partitions=10, sampling="under")
kfold.kfold(rrna_pair_probs, not_rrna_pair_probs, save_folder="/projects/bio/rrna/data/rnafold_results/", n_partitions=10, sampling="over")
