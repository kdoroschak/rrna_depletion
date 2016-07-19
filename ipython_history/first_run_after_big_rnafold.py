# coding: utf-8
get_ipython().magic(u'cd rrna/src')
import numpy as np
mean_pair_probs = np.load("../data/rnafold_results/rnafold_mean_pair_probs.npy")
rrna_pair_probs = np.load("../data/rnafold_results/rnafold_rrna_pair_probs.npy")
mean_pair_probs.shape
rrna_pair_probs.shape
labels = np.zeros(mean_pair_probs.shape[0] + rrna_pair_probs[0])
labels = np.zeros(mean_pair_probs.shape[0] + rrna_pair_probs.shape[0])
labels.shape
for i in range(mean_pair_probs.shape[0], -1):
    print i
    
for i in range(mean_pair_probs.shape[0], labels.shape[0]-1):
    labels[i] = 1
    
np.count_nonzero(labels)
for i in range(mean_pair_probs.shape[0]-1, labels.shape[0]-1):
    labels[i] = 1
    
np.count_nonzero(labels)
pair_probs = np.hstack(mean_pair_probs, rrna_pair_probs)
pair_probs = np.hstack([mean_pair_probs, rrna_pair_probs])
pair_probs = np.vstack([mean_pair_probs, rrna_pair_probs])
import kfold
kfold.kfold(labels, pair_probs)
get_ipython().magic(u'save')
