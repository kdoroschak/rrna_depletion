# coding: utf-8
get_ipython().magic(u'cd rrna/src')
from featureExtraction import *
import generate_rnafold_features as rnafold
fe = FeatureExtractor()
fe.load_chunked_sequences()
results = rnafold.predict_secondary_for_chunks_parallel(fe, "rnafold_mean", gene_list="/homes/gws/kdorosch/unfinished.npy", n_processes=40)
