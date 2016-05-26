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

home_path = os.path.expanduser("~")

def check_chromosome(chrom):
    if chrom in ["chr%d" % (i+1) for i in range(22)]:
        return chrom, True
    elif chrom in ["chrX", "chrY"]:
        return chrom, True
    elif chrom in ["%d" % (i+1) for i in range(22)]:
        return "chr%s" % chrom, True
    elif chrom in ["X", "Y"]:
        return "chr%s" % chrom, True
    else:
        return None, False


class FeatureExtractor(object):
    def __init__(self,
                 fasta_file_path=home_path+"/rrna/data/ref_genomes/GRCh38_o.p3.genome.fa",
                 gtf_file_path=home_path+"/rrna/data/annotations/Homo_sapiens.GRCh38.84.gtf",
                 save_path=home_path+"/rrna/data/featureExtractors/Homo_sapiens.GRCh38.84_gene_features.pkl"):

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
                                        dbfn=self.gtf_db_path,
                                        disable_infer_transcripts=True,
                                        disable_infer_genes=True)

        gtf_features = np.array(list(gtf_db.all_features()))
        gtf_types = np.array([f.featuretype for f in gtf_features])
        # gtf_genes = gtf_features[gtf_types == "gene"]
        gtf_exons = gtf_features[gtf_types == "exon"]
        # gtf_features = None
        # gtf_types = None

        for gtf_exon in gtf_exons:
            chrom, accepted = check_chromosome(gtf_exon.chrom)

            if accepted:
                if "gene_name" in gtf_exon.attributes.keys():
                    gene_name = gtf_exon.attributes["gene_name"][0]
                else:
                    gene_name = gtf_exon.attributes["gene_id"][0]

                if gene_name in self.genefeatures:
                    self.genefeatures[gene_name].exons.append([gtf_exon.start,
                                                               gtf_exon.end])
                else:
                    gf = GeneFeatures(gene_name, chrom)
                    self.genefeatures[gene_name] = gf
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
            print gf_key
            gf = self.genefeatures[gf_key]
            chrom, accepted = check_chromosome(gf.chrom)
            if accepted:
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
    def __init__(self, gene_name, chrom):
        self.chrom = chrom
        self.name = gene_name

        self.exons = []
        self.chunk_size = 0
        self.chunks = []
        self.chunk_locs = []
        self.sequence = ""
        self.coverage = None

    @property
    def start(self):
        self.exons = sorted(self.exons)
        return np.min(self.exons)

    @property
    def end(self):
        self.exons = sorted(self.exons)
        return np.max(self.exons)

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
        self.coverage = np.mean([c.n for c in pileupcolumns])
        return self.coverage

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

    def calculate_kmers(self, k, ignore_N=True):
        kmers = []

        if k < len(self.sequence):
            for pos in range(len(self.sequence) - k):
                seq = self.sequence[pos: pos + k]
                if not "N" in seq:
                    kmers.append(seq)

        return kmers

    def calculate_kmers_and_coverage(self, k, pfile, ignore_N=True):
        kmers_cov = []
        if k < len(self.sequence):
            seq_coverage = np.zeros(self.end - self.start)
            pileupcolumns = pfile.pileup(self.chrom, self.start, self.end)

            for column in pileupcolumns:
                if 0 <= column.pos-self.start < len(seq_coverage):
                    seq_coverage[column.pos-self.start] = column.n

            self.exons = sorted(self.exons)

            glued_seq_coverage = []
            for exon in self.exons:
                glued_seq_coverage += seq_coverage[exon[0]: exon[1]].tolist()

            glued_seq_coverage = np.array(glued_seq_coverage)
            for pos in range(len(glued_seq_coverage)-k):
                seq = self.sequence[pos: pos + k]
                if not ignore_N or not "N" in seq:
                    kmers_cov.append([seq,
                                      np.mean(glued_seq_coverage[pos: pos+k])])

        return kmers_cov