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


def extract_kmers_from_sequence_thread(args):
    seq = args[0]
    k =  args[1]
    rc =  args[2]

    kmers = set()
    for ipos in range(0, len(seq) - k):
        this_seq = seq[ipos: ipos + k]
        kmers = kmers.union([this_seq])
        if rc:
            kmers = kmers.union([this_seq[::-1]])

    return kmers


def extract_kmers_from_sequence(seq, k=6, reverse_complement=False, n_worker=1,
                                work_mulitplicator=3):
    kmers = set()

    if len(seq) < 1e4:
        n_worker = 1

    multi_params = []
    step = int(len(seq) / n_worker / work_mulitplicator)
    for i in range(n_worker * work_mulitplicator):
        if i < n_worker - 1:
            multi_params.append([seq[i * step: (i+1) * step], k, reverse_complement])
        else:
            multi_params.append([seq[i * step: ], k, reverse_complement])

    if n_worker > 1:
        pool = Pool(n_worker)
        results = pool.map(extract_kmers_from_sequence_thread, multi_params)
        pool.close()
        pool.join()
    else:
        results = map(extract_kmers_from_sequence_thread, multi_params)

    for result in results:
        kmers = kmers.union(result)

    return kmers


def extract_kmers_in_ribosome_region(k=6,
                                     ribo_path="/homes/gws/sdorkenw/rrna/data/rrna_hg38.gtf",
                                     reference_fasta_path="/homes/gws/sdorkenw/rrna/data/ref_genomes/GRCh38.p3_genomic.fa",
                                     # reference_chr_path="/homes/gws/sdorkenw/rrna/data/chr_fasta/*.fna",
                                     save_folder="/homes/gws/sdorkenw/rrna/data/",
                                     reload=True, reverse_complement=True,
                                     n_worker=1):
    if os.path.exists(save_folder + "ribo_features.pkl") and reload:
        f = open(save_folder + "ribo_features.pkl", "r")
        ribo_features = pkl.load(f)
        f.close()
    else:
        ribo_db = gffutils.create_db(ribo_path, dbfn=":memory:")
        ribo_features = list(ribo_db.all_features())
        f = open(save_folder + "ribo_features.pkl", "w")
        pkl.dump(ribo_features, f, protocol=pkl.HIGHEST_PROTOCOL)
        f.close()

    if os.path.exists(save_folder + "fasta_seqs.pkl") and reload:
        f = open(save_folder + "fasta_seqs.pkl", "r")
        fasta_seqs = pkl.load(f)
        f.close()
    else:
        # reference_chr_paths = glob.glob(reference_chr_path)
        # fasta_seqs = {}
        # for path in reference_chr_paths:
        #     print(re.findall("[\w\d]+", os.path.basename(path))[0])
        #     fasta_seqs[re.findall("[\w\d]+", os.path.basename(path))[0]] = \
        #         SeqIO.parse(open(path), 'fasta').next()

        fasta_seqs = {}
        fasta_parse = SeqIO.parse(open(reference_fasta_path), 'fasta')
        for anno in fasta_parse:
            fasta_seqs[anno.name] = anno.seq

        f = open(save_folder + "fasta_seqs.pkl", "w")
        pkl.dump(fasta_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)
        f.close()

    cm_ids = ["CM000%d.2" % (663 + i) for i in range(24)]
    cm_mapper = {}
    for icm_id in range(1, 1+len(cm_ids)):
        if icm_id < 23:
            cm_mapper["chr%d" % icm_id] = cm_ids[icm_id-1]
        elif icm_id == 23:
            cm_mapper["chrX"] = cm_ids[icm_id-1]
        else:
            cm_mapper["chrY"] = cm_ids[icm_id-1]

    # print "Number of sequences:", len(ribo_features)
    kmers = set()
    cnt = 0
    unknown_seqs = set()

    for feat in ribo_features:
        # if feat.chrom in fasta_seqs:
        if feat.chrom in cm_mapper:
            chrom = cm_mapper[feat.chrom]
            seq = fasta_seqs[chrom][feat.start: feat.end].tostring()
            if feat.featuretype == "exon":
                # print feat.chrom, "sequence length:", len(seq)
                kmers = kmers.union(extract_kmers_from_sequence(seq, k=k, reverse_complement=reverse_complement,
                                                                n_worker=n_worker))
        else:
            # print "unknown:", feat.chrom
            unknown_seqs = unknown_seqs.union(set([feat.chrom]))
        # print cnt, len(kmers), len(unknown_seqs)
        cnt += 1

    # print len(unknown_seqs)
    return kmers


def kmer_analysis_per_base_thread(args):
    bam_path = args[0]
    ref = args[1]
    ref_length = args[2]
    k = args[3]
    contiguous = args[4]
    reference_seqs = args[5]
    quality_threshold = args[6]

    pfile = pysam.AlignmentFile(bam_path, "rb")
    pileupcolumns = pfile.pileup(ref, 0, ref_length)

    ribo_kmers = extract_kmers_in_ribosome_region(k=k, reverse_complement=True)

    pos = 0
    coverage_ribo = []
    coverage_nonribo = []
    coverage = [0 for _ in range(k)]
    cnt_k = k
    
    time_start = time.time()
    for column in pileupcolumns:
        while pos <= column.pos:
            if cnt_k < k or pos % 1000 == 0 or contiguous:
                if column.n > quality_threshold and pos == column.pos:
                    coverage[pos % k] = column.n
                else:
                    coverage[pos % k] = 0

                this_coverage = np.mean(coverage)
                cnt_k -= 1

            if cnt_k == 0 or contiguous:
                if pos >= k:
                    kmer = reference_seqs[pos - k: pos]
                    if not "N" in kmer:
                        if kmer in ribo_kmers:
                            coverage_ribo.append(this_coverage)
                        else:
                            coverage_nonribo.append(this_coverage)

                cnt_k = k

            if pos % 1e6 == 0 and pos > 0:
                print "%s is at %d of %d; %.3f%%; ~%.2fmin to go" % \
                      (ref, pos, ref_length, float(pos)/ref_length*100,
                       (ref_length / float(pos) - 1) * (time.time()-time_start) / 60)
            pos += 1

    print "%s: finished" % ref
    reference_seqs = None
    pileupcolumns = None
    return coverage_ribo, coverage_nonribo


def kmer_analysis_per_base(bam_path="/homes/gws/sdorkenw/rrna/data/alignments/SRR891244.Aligned.sortedByCoord.out.bam",
                  # fa_path="/homes/gws/sdorkenw/rrna/data/ref_genomes/GRCh38.p3_genomic.fa",
                  fa_path="/homes/gws/sdorkenw/reference_genome_38/GRCh38_o.p3.genome.fa",
                  save_path="/homes/gws/sdorkenw/rrna/data/kmer_analysis/",
                  quality_threshold=0, contiguous=False, k=14, n_workers=10):
    pfile = pysam.AlignmentFile(bam_path, "rb")
    references = []
    ref_lengths = {}

    for ir in range(len(pfile.references)):
        ref = pfile.references[ir]
        if ref.startswith("chr"):
            references.append(ref)
            ref_lengths[ref] = pfile.lengths[ir]

    print "create reference sequences"
    reference_seqs = {}
    for rec in SeqIO.parse(open(fa_path), 'fasta'):
        if rec.id.startswith("chr"):
            reference_seqs[rec.id] = rec.seq

    multi_params = []
    for ref in references:
        multi_params.append([bam_path, ref, ref_lengths[ref], k, contiguous,
                             reference_seqs[ref], quality_threshold])


    if n_workers > 1:
        pool = Pool(n_workers)
        results = pool.map(kmer_analysis_per_base_thread, multi_params)
        pool.close()
        pool.join()
    else:
        results = map(kmer_analysis_per_base_thread, multi_params)

    coverage_nonribo = []
    coverage_ribo = []
    for result in results:
        coverage_ribo.append(result[0])
        coverage_nonribo.append(result[1])

    if contiguous:
        np.save(save_path + "coverage_ribo_k%d_%s_cont" % (k, re.findall("[\d]+", bam_path)[-2]), coverage_ribo)
        np.save(save_path + "coverage_nonribo_k%d_%s_cont" % (k, re.findall("[\d]+", bam_path)[-2]), coverage_nonribo)
    else:
        np.save(save_path + "coverage_ribo_k%d_%s" % (k, re.findall("[\d]+", bam_path)[-2]), coverage_ribo)
        np.save(save_path + "coverage_nonribo_k%d_%s" % (k, re.findall("[\d]+", bam_path)[-2]), coverage_nonribo)

def kmer_analysis_per_read_thread(args):
    reads = args[0]
    k = args[1]
    start_only = args[2]

    ribo_kmers = extract_kmers_in_ribosome_region(k=k, reverse_complement=True)

    pos = 0
    reads_ribo = {}
    reads_nonribo = {}
    cnt = 0

    time_start = time.time()
    for read in reads:
        seq = read.seq
        if start_only:
            seq_length = k+1
        else:
            seq_length = len(seq)

        for pos in range(0, seq_length-k):
            kstring = ''.join(seq[pos: pos+k])
            if seq[pos: pos+k] in ribo_kmers:
                if kstring in reads_ribo:
                    reads_ribo[kstring] += 1
                else:
                    reads_ribo[kstring] = 1
            else:
                if kstring in reads_nonribo:
                    reads_nonribo[kstring] += 1
                else:
                    reads_nonribo[kstring] = 1

        cnt += 1
        if cnt % 5e4 == 0 and cnt > 0:
            print "I am at %d of %d; %.3f%%; ~%.2fmin to go" % \
                  (cnt, len(reads), float(cnt)/len(reads)*100,
                   (len(reads) / float(cnt) - 1) * (time.time()-time_start) / 60)

    print "finished"
    reads = None
    ribo_kmers = None
    return reads_nonribo, reads_ribo


def kmer_analysis_per_read(bam_path="/homes/gws/sdorkenw/rrna/data/alignments/SRR891244.Aligned.sortedByCoord.out.bam",
                  save_path="/homes/gws/sdorkenw/rrna/data/kmer_analysis/",
                  start_only=False, k=14, n_workers=10):
    pfile = pysam.AlignmentFile(bam_path, "rb")

    split_every = (lambda n, it: takewhile(bool, (list(islice(it, n)) for _ in count(0))))
    split_reads = split_every(np.ceil(pfile.mapped/float(n_workers*3)), pfile.fetch())

    multi_params = []
    for _ in range(n_workers):
        multi_params.append([split_reads.next(), k, start_only])

    if n_workers > 1:
        pool = Pool(n_workers)
        results = pool.map(kmer_analysis_per_read_thread, multi_params)
        pool.close()
        pool.join()
    else:
        results = map(kmer_analysis_per_read_thread, multi_params)

    kmers_nonribo = {}
    kmers_ribo = {}
    for result in results:
        for key in result[0].keys():
            if key in kmers_nonribo:
                kmers_nonribo[key] += result[0][key]
            else:
                kmers_nonribo[key] = result[0][key]
        for key in result[1].keys():
            if key in kmers_ribo:
                kmers_ribo[key] += result[1][key]
            else:
                kmers_ribo[key] = result[1][key]

    if start_only:
        np.save(save_path + "read_ribo_k%d_%s_start_only" % (k, re.findall("[\d]+", bam_path)[-1]), kmers_ribo)
        np.save(save_path + "read_nonribo_k%d_%s_start_only" % (k, re.findall("[\d]+", bam_path)[-1]), kmers_nonribo)
    else:
        np.save(save_path + "read_ribo_k%d_%s" % (k, re.findall("[\d]+", bam_path)[-1]), kmers_ribo)
        np.save(save_path + "read_nonribo_k%d_%s" % (k, re.findall("[\d]+", bam_path)[-1]), kmers_nonribo)