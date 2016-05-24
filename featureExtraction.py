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
import scipy.spatial
import time


def generate_chrom_mapper(factor=1e6):
    chrom_mapper = {}

    rev_chrom_mapper = {}
    for ichrom in range(22):
        chrom_mapper["chr%d" % (ichrom + 1)] = ichrom * factor
        rev_chrom_mapper[ichrom * factor] = "chr%d" % (ichrom + 1)

    chrom_mapper["chrX"] = 23 * factor
    rev_chrom_mapper[23 * factor] = "chrX"
    chrom_mapper["chrY"] = 24 * factor
    rev_chrom_mapper[24 * factor] = "chrY"

    return chrom_mapper


class FeatureExtractor(object):
    def __init__(self, fasta_file_path, gtf_file_path, save_path):
        self.fasta_file_path = fasta_file_path
        self.gtf_file_path = gtf_file_path
        self.gtf_db_path = self.gtf_file_path[:-4] + "features.db"
        self.genefeatures = {}
        self.gf_save_path = save_path

    def load_gene_features(self):
        if os.path.exists(self.gf_save_path):
            try:
                f = open(self.gf_save_path, "r")
                self.genefeatures = pkl.load(f)
                f.close()
            except:
                return False
        else:
            return False
        return True

    def generate_gene_features(self, reload=True):
        if os.path.exists(self.gtf_db_path) and reload:
            gtf_db = gffutils.FeatureDB(self.gtf_db_path)
        else:
            if os.path.exists(self.gtf_db_path):
                os.remove(self.gtf_db_path)
            gtf_db = gffutils.create_db(self.gtf_file_path,
                                        dbfn=self.gtf_db_path)

        gtf_features = np.array(list(gtf_db.all_features()))
        gtf_types = np.array([f.featuretype for f in gtf_features])
        gtf_genes = gtf_features[gtf_types == "gene"]
        gtf_exons = gtf_features[gtf_types == "exon"]
        # gtf_features = None
        # gtf_types = None

        chrom_mapper = generate_chrom_mapper()

        gene_coordinates = []
        name_list = []
        for gtf_gene in gtf_genes:
            if gtf_gene.chrom in chrom_mapper:
                gf = GeneFeatures(gtf_gene)
                self.genefeatures[gtf_gene.id] = gf
                gene_coordinates.append([chrom_mapper[gf.chrom], gf.start])
                gene_coordinates.append([chrom_mapper[gf.chrom], gf.end])
                name_list.append(gf.name)
                name_list.append(gf.name)

        raise()
        gene_kdtree = scipy.spatial.cKDTree(gene_coordinates)

        for gtf_exon in gtf_exons:
            if gtf_exon.chrom in chrom_mapper:
                _, id = gene_kdtree.query([chrom_mapper[gtf_exon.chrom],
                                           gtf_exon.start])
                gene_name = name_list[id]
                assert gene_coordinates[id][0] == \
                       chrom_mapper[self.genefeatures[gene_name].chrom]
                self.genefeatures[gene_name].exons.append([gtf_exon.start,
                                                           gtf_exon.end])

        f = open(self.gf_save_path, "w")
        pkl.dump(self.genefeatures, f)
        f.close()

    def generate_chunked_sequences(self):
        fasta_seqs = {}
        fasta_parse = SeqIO.parse(open(self.fasta_file_path), 'fasta')
        for anno in fasta_parse:
            fasta_seqs[anno.name] = anno.seq

        for gf_key in self.genefeatures.keys():
            gf = self.genefeatures[gf_key]
            gf.generate_full_sequence(fasta_seqs[gf.chrom])
            gf.partition_sequence_in_chunks(chunk_size=50)

        # f = open(self.gf_save_path, "w")
        # pkl.dump(self.genefeatures, f)
        # f.close()

    def generate_sequences(self):
        fasta_seqs = {}
        fasta_parse = SeqIO.parse(open(self.fasta_file_path), 'fasta')
        for anno in fasta_parse:
            fasta_seqs[anno.name] = anno.seq

        for gf_key in self.genefeatures.keys():
            gf = self.genefeatures[gf_key]
            gf.generate_full_sequence(fasta_seqs[gf.chrom])

    def calculate_gene_coverage(self, bam_path):
        pfile = pysam.AlignmentFile(bam_path, "rb")

        coverage_dict = {}
        for gf_key in self.genefeatures.keys():
            gf = self.genefeatures[gf_key]
            coverage_dict[gf_key] = gf.calculate_chunkwise_coverage(pfile)

        return coverage_dict


class GeneFeatures(object):
    def __init__(self, gtf_gene):
        self.chrom = gtf_gene.chrom
        self.start = gtf_gene.start
        self.end = gtf_gene.end
        self.name = gtf_gene.id

        self.exons = []
        self.chunk_size = 0
        self.chunks = []
        self.chunk_locs = []
        self.sequence = ""

    def generate_full_sequence(self, fasta_seq):
        self.exons = sorted(self.exons)
        self.sequence = ""
        for exon in self.exons:
            self.sequence += fasta_seq[exon[0]: exon[1]].tostring()

    def partition_sequence_in_chunks(self, chunk_size=50):
        assert chunk_size > 0
        assert chunk_size % 2 == 0

        self.chunk_size = chunk_size

        if len(self.sequence) >= chunk_size:
            self.chunks = []
            pos = 0
            while pos < len(self.sequence):
                if len(self.sequence)-pos < chunk_size:
                    self.chunks.append(self.sequence[-chunk_size:])
                    self.chunk_locs.append([len(self.sequence) - chunk_size,
                                            len(self.sequence)])
                    pos = len(self.sequence)
                else:
                    self.chunks.append(self.sequence[pos: pos+chunk_size])
                    self.chunk_locs.append([pos, pos+chunk_size])
                    pos += chunk_size/2

    def calculate_mean_coverage(self, pfile):
        pileupcolumns = pfile.pileup(self.chrom, self.start, self.end)
        return np.mean([c.n for c in pileupcolumns])

    def calculate_chunkwise_coverage(self, pfile):
        coverage = []
        pileupcolumns = pfile.pileup(self.chrom, self.start, self.end)
        gene_coverage = {}
        for column in pileupcolumns:
            gene_coverage[column.pos] = column.n

        for chunk_loc in self.chunk_locs:
            chunk_coverage = 0
            for pos in range(chunk_loc[0], chunk_loc[1]):
                if pos in gene_coverage:
                    chunk_coverage += gene_coverage[pos]

            coverage.append(chunk_coverage / float(self.chunk_size))

        return coverage

    def calculate_kmers(self, k):
        kmers = []

        if k > len(self.sequence):
            for pos in range(len(self.sequence) - k):
                kmers.append(self.sequence[pos:, pos + k])

        return kmers

    def calculate_kmers_and_coverage(self, k, pfile):
        kmers_cov = []
        if k > len(self.sequence):
            seq_coverage = np.zeros(len(self.sequence))
            pileupcolumns = pfile.pileup(self.chrom, self.start, self.end)

            for column in pileupcolumns:
                seq_coverage[column.pos-self.start] = column.n

            for pos in range(len(self.sequence)-k):
                kmers_cov.append([self.sequence[pos:, pos+k],
                                  np.mean(seq_coverage[pos: pos+k])])

        return kmers_cov