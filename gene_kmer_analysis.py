import gffutils # pip install gffutils
from Bio import SeqIO # pip install Biopython
import pysam # pip install pysam
import cPickle as pkl
from itertools import islice, takewhile, count
import glob
import multiprocessing
import multiprocessing.pool
import numpy as np
import os
import plotUtils as pu
import re
import time

import featureExtraction


def extract_ribo_kmers(k=18,
                       gtf_path="/homes/gws/sdorkenw/rrna/data/ref_genomes/rrna_hg38.gtf",
                       fasta_file_path="/homes/gws/sdorkenw/reference_genome_38/GRCh38_o.p3.genome.fa",
                       save_path="/homes/gws/sdorkenw/rrna/data/featureExtractors/ribo_features.db"):
    fe = featureExtraction.FeatureExtractor(fasta_file_path, gtf_path,
                                            save_path)

    if not fe.load_gene_features():
        fe.generate_gene_features()

    fe.generate_chunked_sequences()

    kmers = set()
    for gf_key in fe.genefeatures.keys():
        gf = fe.genefeatures[gf_key]
        kmers = kmers.union(set(gf.calculate_kmers(k=18)))

    return kmers





