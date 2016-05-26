import gffutils # pip install gffutils
from Bio import SeqIO # pip install Biopython
import pysam # pip install pysam
import cPickle as pkl
from itertools import islice, takewhile, count
import glob
from multiprocessing import Pool
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
    # fe = featureExtraction.FeatureExtractor(fasta_file_path, gtf_path,
    #                                         save_path)
    fe = featureExtraction.FeatureExtractor()

    if not fe.load_gene_features():
        fe.generate_gene_features(reload=True)

    # fe.generate_chunked_sequences()

    kmers = set()
    for gf_key in fe.genefeatures.keys():
        print gf_key
        gf = fe.genefeatures[gf_key]
        kmers = kmers.union(set(gf.calculate_kmers(k=k)))

    return kmers


def _extract_gene_kmers_thread(args):
    bam_path = args[0]
    k = args[1]
    gf = args[2]

    pfile = pysam.AlignmentFile(bam_path, "rb")

    return gf.name, gf.calculate_kmers_and_coverage(k=k, pfile=pfile)


def extract_gene_kmers(bam_path,
                       k=18,
                       n_worker=30,
                       gtf_path="/homes/gws/sdorkenw/rrna/data/annotations/Homo_sapiens.GRCh38.84.gtf",
                       fasta_file_path="/homes/gws/sdorkenw/reference_genome_38/GRCh38_o.p3.genome.fa",
                       save_path="/homes/gws/sdorkenw/rrna/data/featureExtractors/gene_features.pkl"):
    # fe = featureExtraction.FeatureExtractor(fasta_file_path, gtf_path,
    #                                         save_path)
    fe = featureExtraction.FeatureExtractor()

    if not fe.load_gene_features():
        fe.generate_gene_features(reload=True)

    fe.generate_chunked_sequences()

    # pfile = pysam.AlignmentFile(bam_path, "rb")
    kmers_cov = {}

    multi_params = []
    for gf_key in fe.genefeatures.keys():
        multi_params.append([bam_path, k, fe.genefeatures[gf_key]])

    if n_worker > 1:
        pool = Pool(n_worker)
        results = pool.map_async(_extract_gene_kmers_thread, multi_params)
        pool.close()  # No more work
        while True:
            if results.ready():
                break
            remaining = results._number_left
            print "Waiting for", remaining, "tasks to complete..."
            time.sleep(1.)
        pool.join()
        results = results.get()
    else:
        results = map(_extract_gene_kmers_thread, multi_params)

    for result in results:
        kmers_cov[result[0]] = result[1]


    return kmers_cov


def kmer_analysis_gene_level(k, bam_path, reload=True,
                             ribo_gene_name_path="/homes/gws/sdorkenw/rrna/data/annotations/Homo_sapiens.GRCh38.84_rrna_gene_names.npy",
                             gtf_path="/homes/gws/sdorkenw/rrna/data/annotations/Homo_sapiens.GRCh38.84.gtf",
                             fasta_file_path="/homes/gws/sdorkenw/reference_genome_38/GRCh38_o.p3.genome.fa",
                             save_path="/homes/gws/sdorkenw/rrna/data/featureExtractors/gene_features.pkl"):

    kmer_save_path = "/homes/gws/sdorkenw/rrna/data/featureExtractors/kmer_cov_%s_%d.pkl" % (os.path.basename(bam_path).split("_")[0], k)
    if os.path.exists(kmer_save_path) and reload:
        f = open(kmer_save_path, "r")
        kmers_cov = pkl.load(f)
        f.close()
    else:
        kmers_cov = extract_gene_kmers(bam_path, k, n_worker=50)

        f = open(kmer_save_path, "w")
        pkl.dump(kmers_cov, f)
        f.close()

    print "Got kmers"
    ribo_names = np.load(ribo_gene_name_path)
    nonribo_names = []

    ribo_kmers = set()

    cnt = 0
    kmers_cov_keys = kmers_cov.keys()
    for igf_key in range(len(kmers_cov_keys)):
        print "%d of %d" % (igf_key, len(kmers_cov_keys))
        gf_key = kmers_cov_keys[igf_key]
        if gf_key in ribo_names:
            cnt += 1
            if len(kmers_cov[gf_key]) > 0:
                ribo_kmers = ribo_kmers.union(kmers_cov[gf_key][:, 0])
        else:
            nonribo_names.append(gf_key)

    nonzero = []
    for key in kmers_cov.keys():
        if len(kmers_cov[key]) > 0:
            nonzero.append(kmers_cov[key])
    raise()
    ribo_bin = []
    nonribo_bin = []

    for inonribo_name in range(len(nonribo_names)):
        print "%d of %d" % (inonribo_name, len(nonribo_names))
        nonribo_name = nonribo_names[inonribo_name]
        for kmer in kmers_cov[nonribo_name]:
            cov = int(kmer[1])
            if kmer[0] in ribo_kmers:
                ribo_bin.append(cov)
            else:
                nonribo_bin.append(cov)

    return ribo_bin, nonribo_bin



# def kmer_analysis_gene_level