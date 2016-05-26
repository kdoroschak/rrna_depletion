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
                       save_path="/homes/gws/sdorkenw/rrna/data/featureExtractors/ribo_features.pkl"):
    fe = featureExtraction.FeatureExtractor(fasta_file_path, gtf_path,
                                            save_path)

    if not fe.load_gene_features():
        fe.generate_gene_features(reload=True)

    # fe.generate_chunked_sequences()

    kmers = set()
    for gf_key in fe.genefeatures.keys():
        gf = fe.genefeatures[gf_key]
        kmers = kmers.union(set(gf.calculate_kmers(k=k)))

    return kmers


def extract_gene_kmers(bam_path, k=18,
                       gtf_path="/homes/gws/sdorkenw/rrna/data/annotations/Homo_sapiens.GRCh38.84.gtf",
                       fasta_file_path="/homes/gws/sdorkenw/reference_genome_38/GRCh38_o.p3.genome.fa",
                       save_path="/homes/gws/sdorkenw/rrna/data/featureExtractors/gene_features.pkl"):
    fe = featureExtraction.FeatureExtractor(fasta_file_path, gtf_path,
                                            save_path)

    # if not fe.load_gene_features():
    fe.generate_gene_features(reload=True)

    # fe.generate_chunked_sequences()

    pfile = pysam.AlignmentFile(bam_path, "rb")
    kmers_cov = {}
    for gf_key in fe.genefeatures.keys():
        print gf_key
        gf = fe.genefeatures[gf_key]
        kmers_cov[gf_key] = np.array(gf.calculate_kmers_and_coverage(k=k, pfile=pfile))

    return kmers_cov


def kmer_analysis_gene_level(k, bam_path,
                             ribo_gene_name_path="/homes/gws/sdorkenw/rrna/data/annotations/Homo_sapiens.GRCh38.84_rrna_gene_names.gtf",
                             gtf_path="/homes/gws/sdorkenw/rrna/data/annotations/Homo_sapiens.GRCh38.84.gtf",
                             fasta_file_path="/homes/gws/sdorkenw/reference_genome_38/GRCh38_o.p3.genome.fa",
                             save_path="/homes/gws/sdorkenw/rrna/data/featureExtractors/gene_features.pkl"):

    kmers_cov = extract_gene_kmers(bam_path, k, gtf_path, fasta_file_path,
                                   save_path)

    f = open("/homes/gws/sdorkenw/rrna/data/featureExtractors/kmer_cov_SRR891241_18.pkl", "w")
    pkl.dump(kmers_cov, f)
    f.close()

    print "Got kmers"
    ribo_names = np.load(ribo_gene_name_path)
    nonribo_names = []

    ribo_kmers = set()

    for gf_key in kmers_cov.keys():
        if gf_key in ribo_names:
            ribo_kmers = ribo_kmers.union(kmers_cov[gf_key][:, 0])
        else:
            nonribo_names.append(gf_key)

    ribo_bin = []
    nonribo_bin = []

    for nonribo_name in nonribo_names:
        for kmer in kmers_cov[nonribo_name]:
            if kmer[0] in ribo_kmers:
                ribo_bin.append(kmer[1])
            else:
                nonribo_bin.append(kmer[1])

    return ribo_bin, nonribo_bin


