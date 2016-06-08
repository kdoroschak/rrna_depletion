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
import scipy.spatial
import sklearn.neighbors
import sklearn.decomposition
import sklearn.preprocessing
import sklearn.ensemble
import sklearn.linear_model
import time
import rna_analysis as ra

import featureExtraction

home_path = os.path.expanduser("~")

def extract_ribo_kmers(k=18,
                       gtf_path="/homes/gws/sdorkenw/rrna/data/ref_genomes/rrna_hg38.gtf",
                       fasta_file_path="/homes/gws/sdorkenw/reference_genome_38/GRCh38_o.p3.genome.fa",
                       save_path="/homes/gws/sdorkenw/rrna/data/featureExtractors/ribo_features.pkl"):
    # fe = featureExtraction.FeatureExtractor(fasta_file_path, gtf_path,
    #                                         save_path)
    fe = featureExtraction.FeatureExtractor()

    if not fe.load_gene_features():
        fe.generate_gene_features(reload=True)

    fe.generate_chunked_sequences()

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
                       n_worker=10,
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
            time.sleep(10)
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
        kmers_cov = extract_gene_kmers(bam_path, k, n_worker=10)

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
                ribo_kmers = ribo_kmers.union([k[0] for k in kmers_cov[gf_key]])
        else:
            nonribo_names.append(gf_key)

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

    print "ribo: %.3f with %d elements" % (np.mean(ribo_bin), len(ribo_bin))
    print "nonribo: %.3f with %d elements" % (np.mean(nonribo_bin), len(nonribo_bin))

    return ribo_bin, nonribo_bin


def extract_kmers_from_ribo_genes(k=50, ribo_gene_name_path="/homes/gws/sdorkenw/rrna/data/annotations/Homo_sapiens.GRCh38.84_rrna_gene_names.npy",
                                  as_integers=True):
    ribo_names = np.load(ribo_gene_name_path)
    ribo_kmers = set()

    fe = featureExtraction.FeatureExtractor()

    if not fe.load_gene_features():
        fe.generate_gene_features(reload=True)

    fe.generate_chunked_sequences()

    for gf_key in ribo_names:
        print gf_key
        if gf_key in fe.genefeatures:
            gf = fe.genefeatures[gf_key]
            ribo_kmers = ribo_kmers.union(set(gf.calculate_kmers(k=k)))

    np.save(home_path + "/rrna/data/featureExtractors/ribo_kmers_from_genes_%d.npy" % k,
            list(ribo_kmers))

    if as_integers:
        mapper = {"A": 1, "C": 2, "G": 3, "T": 4}
        ribo_kmers_integers = np.zeros([len(ribo_kmers), k])
        i_kmer = 0
        for kmer in ribo_kmers:
            print "%d of %d" % (i_kmer, len(ribo_kmers))
            ribo_kmers_integers[i_kmer] = np.array([mapper[nuc] for nuc in kmer], dtype=np.int)
            i_kmer += 1

        np.save(home_path + "/rrna/data/featureExtractors/ribo_kmers_from_genes_%d_int.npy" % k,
                ribo_kmers_integers)

def extract_expression_features(bam_path,
                                path_ordered_gene_names = home_path +
                                "/rrna/data/featureExtractors/Homo_sapiens.GRCh38.84_ordered_gene_names.npy",
                                n_worker=10):
    gene_names = np.load(path_ordered_gene_names)

    identifier = re.findall("[\d]+", os.path.basename(bam_path))[0]
    multi_params = []
    for i_worker in range(n_worker):
        save_path = home_path + "/rrna/data/featureExtractors/coverage_features_%s_%d.pkl" % (identifier, i_worker)
        multi_params.append([gene_names[i_worker::n_worker], bam_path, save_path])

    if n_worker > 1:
        pool = Pool(n_worker)
        pool.map(_extract_expression_features_thread, multi_params)
        pool.close()
        pool.join()
    else:
        map(_extract_expression_features_thread, multi_params)


def _extract_expression_features_thread(args):
    gene_names = args[0]
    bam_path = args[1]
    save_path = args[2]

    fe = featureExtraction.FeatureExtractor()

    if not fe.load_gene_features():
        fe.generate_gene_features(reload=True)

    fe.generate_chunked_sequences()

    if os.path.exists(save_path):
        with open(save_path, "r") as f:
            coverage = pkl.load(f)
    else:
        coverage = {}

    pfile = pysam.AlignmentFile(bam_path, "rb")

    time_start = time.time()
    for igf_key in range(len(gene_names)):
        print igf_key, time.time()-time_start
        gf_key = gene_names[igf_key]
        if not gf_key in coverage and gf_key in fe.genefeatures:
            if len(fe.genefeatures[gf_key].chunks) > 0:
                fe.genefeatures[gf_key].calculate_chunkwise_coverage(pfile)
                # print fe.genefeatures[gf_key].chunk_coverage
                try:
                    mean = np.mean(fe.genefeatures[gf_key].chunk_coverage, axis=1)
                    std = np.std(fe.genefeatures[gf_key].chunk_coverage, axis=1)
                    gap = np.max(fe.genefeatures[gf_key].chunk_coverage, axis=1) - \
                          np.min(fe.genefeatures[gf_key].chunk_coverage, axis=1)
                except IndexError:
                    mean = np.zeros(len(fe.genefeatures[gf_key].chunks))
                    std = np.zeros(len(fe.genefeatures[gf_key].chunks))
                    gap = np.zeros(len(fe.genefeatures[gf_key].chunks))
                coverage[gf_key] = np.vstack([mean.T, std.T, gap.T]).T

            if igf_key % 10 == 0:
                with open(save_path, "w") as f:
                    pkl.dump(coverage, f)


def combine_expression_features(re_path):
    paths = glob.glob(re_path)
    exp_dict = {}
    for path in paths:
        print path
        with open(path) as f:
            this_dict = pkl.load(f)

        exp_dict.update(this_dict)

    save_path = paths[0][:-6] + ".pkl"
    print save_path
    with open(save_path, "w") as f:
        pkl.dump(exp_dict, f)

    # return exp_dict



def extract_hamming_distance_feature(path_ribo_kmers=home_path +
                                     "/rrna/data/featureExtractors/ribo_kmers_from_genes_50_int.npy",
                                     path_ordered_gene_names = home_path +
                                     "/rrna/data/featureExtractors/Homo_sapiens.GRCh38.84_ordered_gene_names.npy"):
    fe = featureExtraction.FeatureExtractor()

    if not fe.load_gene_features():
        fe.generate_gene_features(reload=True)

    fe.generate_chunked_sequences()

    nn = sklearn.neighbors.NearestNeighbors(n_neighbors=1, metric="hamming",
                                               n_jobs=30)
    nn.fit(np.load(path_ribo_kmers))

    gene_names = np.load(path_ordered_gene_names)

    save_path = home_path + "/rrna/data/featureExtractors/distance_features.pkl"
    if os.path.exists(save_path):
        with open(save_path, "r") as f:
            distances = pkl.load(f)
    else:
        distances = {}
    mapper = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 5}

    time_start = time.time()
    for igf_key in range(len(gene_names)):
        gf_key = gene_names[igf_key]
        print "%d of %d, time: %.3f" % \
              (igf_key, len(gene_names), time.time() - time_start)

        if not gf_key in distances and gf_key in fe.genefeatures:
            print "next length: %d" % (len(fe.genefeatures[gf_key].chunks))
            chunks = [np.array([mapper[nuc] for nuc in seq], dtype=np.int)
                      for seq in fe.genefeatures[gf_key].chunks]
            distances[gf_key] = nn.kneighbors(chunks)[0]

            if igf_key % 10 == 0:
                with open(save_path, "w") as f:
                    pkl.dump(distances, f)


def agglomerate_features(folder=home_path + "/rrna/data/featureExtractors/",
                         ribo_gene_name_path="/homes/gws/sdorkenw/rrna/data/annotations/Homo_sapiens.GRCh38.84_rrna_gene_names.npy"):
    ribo_genes = np.load(ribo_gene_name_path)
    consensus_genes = set()
    with open(folder + "/distance_features.pkl") as f:
        dist_dict = pkl.load(f)

    consensus_genes = consensus_genes.union(dist_dict.keys())

    expressions = []
    for i in range(1, 7):
        with open(folder + "/coverage_features_89124%d.pkl" % i, "r") as f:
            print folder + "/coverage_features_89124%d.pkl" % i
            expressions.append(pkl.load(f))

        consensus_genes = consensus_genes.intersection(expressions[-1].keys())

    raw_names = np.load(home_path + "/rrna/data/rnafold_results/row_names.npy")[:, 0]
    folding_energy = np.load(home_path + "/rrna/data/rnafold_results/mfe.npy")
    contact_probs = np.load(home_path + "/rrna/data/rnafold_results/pair_probs.npy")

    # folding_dict = {}
    # for i_sample in range(len(folding_energy)):
    #     if not raw_names[i_sample] in folding_dict:
    #         folding_dict[raw_names[i_sample]] = []
    #     folding_dict[raw_names[i_sample]].append([folding_energy[i_sample]] + contact_probs[i_sample].tolist())

    # consensus_genes = consensus_genes.intersection(folding_dict.keys())

    chunks = []
    ribo_chunks = []
    ribo = []
    gene_dists = []
    gene_names = []
    ribo_chunkwise_gene_names = []
    chunkwise_gene_names = []
    gene_cov = []
    print "There are %d genes in the consensus set" % len(consensus_genes)
    for gene in consensus_genes:
        # expressions = expression_dict[gene]
        this_gene_cov = [np.mean(expressions[i][gene][:, 0]) for i in range(6)]

        if np.sum(this_gene_cov) > 0:
            # foldings = folding_dict[gene]
            dists = dist_dict[gene]
            # if len(dists) == len(foldings):
            if True:
                gene_cov.append(this_gene_cov)
                gene_names.append(gene)

                gene_dists.append(np.mean(dists))

                gene_expressions = np.hstack([expressions[i][gene] for i in range(6)])
                relative_expression = []
                for i in range(6):
                    if this_gene_cov[i] > 0:
                        relative_expression.append(expressions[i][gene][:, 0] / this_gene_cov[i])
                    else:
                        relative_expression.append(np.zeros_like(expressions[i][gene][:, 0]))
                relative_expression = np.array(relative_expression).T

                if gene in ribo_genes or gene.startswith("MRP") or gene.startswith("RP"):
                    ribo.append(True)
                    # ribo_chunks += np.hstack([dists, foldings, gene_expressions, relative_expression, np.array([1]*len(foldings)).reshape(len(foldings), 1)]).tolist()
                    ribo_chunks += np.hstack([dists, gene_expressions, relative_expression]).tolist()
                    ribo_chunkwise_gene_names += [gene] * len(dists)
                else:
                    ribo.append(False)
                    # chunks += np.hstack([dists, foldings, gene_expressions, relative_expression, np.array([0]*len(foldings)).reshape(len(foldings))]).tolist()
                    chunks += np.hstack([dists, gene_expressions, relative_expression]).tolist()
                    chunkwise_gene_names += [gene] * len(dists)
            else:
                print gene, len(dists)
        else:
            print gene

    output = [np.array(chunks), np.array(ribo_chunks), np.array(ribo),
              np.array(gene_dists), np.array(gene_cov), np.array(gene_names),
              np.array(chunkwise_gene_names), np.array(ribo_chunkwise_gene_names)]
    np.save(folder + "/agglomerated_features.pkl", output)

    return np.array(chunks), np.array(ribo_chunks), np.array(ribo), \
           np.array(gene_dists), np.array(gene_cov), np.array(gene_names), \
           np.array(chunkwise_gene_names), np.array(ribo_chunkwise_gene_names)



def pca_analysis(data):
    data = sklearn.preprocessing.scale(data)
    pca = sklearn.decomposition.PCA(n_components=.999)
    return pca.fit_transform(data)


def calc_ribo_ratio(count_path, ribo_gene_name_path="/homes/gws/sdorkenw/rrna/data/annotations/Homo_sapiens.GRCh38.84_rrna_gene_names.npy"):
    ribo_genes = np.load(ribo_gene_name_path)

    count_dict = ra.read_values_from_file(count_path)

    with open(home_path + "/rrna/data/featureExtractors/gene_name_mapper.pkl", "r") as f:
        mapper = pkl.load(f)

    non_ribo_count_sum = 0
    ribo_count_sum = 0
    for name in count_dict:
        if name in mapper:
            gene = mapper[name][0]
            if gene in ribo_genes or gene.startswith("MRP") or gene.startswith("RP"):
                ribo_count_sum += count_dict[name]
            else:
                non_ribo_count_sum += count_dict[name]
            # print ribo_count_sum, non_ribo_count_sum

    print "Ribosomal ratio: %.4f" % (ribo_count_sum / (ribo_count_sum + non_ribo_count_sum))


def balance_data(x, y, nbins=10):
    bins = np.linspace(0.0, 0.5, nbins)
    idx = np.digitize(y, bins)
    sums = np.array([np.sum(np.ones_like(y)[idx == k+1]) for k in range(nbins)])
    sums = 1./sums
    sums /= np.min(sums)
    sums[0] = 0
    weights = np.array([sums[idx[y_sample]-1] for y_sample in range(len(y))], dtype=np.int)

    new_x = []
    new_y = []

    for i_sample in range(len(y)):
        new_x += [x[i_sample]]*weights[i_sample]
        new_y += [y[i_sample]]*weights[i_sample]

    return new_x, new_y


def balance_data_binary(x, y):
    y_1 = y[y==1]
    x_1 = x[y==1]
    new_y = y.tolist()
    new_x = x.tolist()
    for i_iter in range(int(np.sum(y == 0)-np.sum(y==0))):
        new_y.append(y_1[i_iter % len(y_1)])
        new_x.append(x_1[i_iter % len(y_1)])
    return new_x, new_y


def regression(x, y, classifier="rfc"):
    # bx, by = balance_data(x, y)

    ids = np.arange(len(y))
    np.random.shuffle(ids)
    bx_tr, by_tr = balance_data_binary(x[ids][:int(.75 * len(x))],
                                       y[ids][:int(.75 * len(x))])
    bx_va, by_va = balance_data_binary(x[ids][:int(.75 * len(x))],
                                       y[ids][:int(.75 * len(x))])

    if classifier == "rfc":
        clf = sklearn.ensemble.RandomForestRegressor(n_jobs=50, n_estimators=20)
    elif classifier == "gbr":
        clf = sklearn.ensemble.GradientBoostingRegressor()
    elif "lin":
        clf = sklearn.linear_model.LinearRegression(normalize=True, n_jobs=50)
    else:
        raise()

    fit = clf.fit(bx_tr, by_tr)

    with open(home_path + "/rrna/data/regression/%s.pkl" % classifier, "w") as f:
        pkl.dump(clf, f)

    pred_y = clf.predict(bx_va)
    print np.mean(np.abs(by_va - pred_y)), np.mean(np.abs(by_va - np.mean(by_va)))
    print np.sqrt(np.mean((pred_y - by_va)**2)), np.sqrt(np.mean((np.mean(by_va) - by_va)**2))
    pred_y = clf.predict(bx)
    print np.mean(np.abs(by-pred_y)), np.mean(np.abs(by-np.mean(by)))
    print np.sqrt(np.mean((pred_y - by)**2)), np.sqrt(np.mean((np.mean(by) - by)**2))


def secondary_structure_distance(non_ribo_structures, ribo_structures):
    kdtree = scipy.spatial.cKDTree(ribo_structures)

    dists = []
    time_start = time.time()
    for cnt in range(len(non_ribo_structures)):
        print "%d - Time: %.3fs" % (cnt, time.time() - time_start)
        cp = kdtree.query(non_ribo_structures[cnt], k=1, n_jobs=50)
        dists.append(cp[0])

    return dists

