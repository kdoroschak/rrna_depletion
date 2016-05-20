import plotUtils as pu
import gffutils
import glob
import numpy as np
import os
import pysam
import re
import scipy.stats


def read_values_from_file(path):
    gene_dict = {}

    with open(path) as f:
        for line in f.readlines():
            gene_name, value = line.split("\t")
            # print gene_name, value
            gene_dict[gene_name] = float(value)

    return gene_dict


def rank_genes(gene_dict, invert_sorting=True):
    keys = gene_dict.keys()
    values = gene_dict.values()

    order = np.argsort(values)
    ranks = np.empty(len(order), np.int)
    if invert_sorting:
        ranks[order] = np.arange(len(order))[::-1]
    else:
        ranks[order] = np.arange(len(order))
    ranking_dict = dict(zip(keys, ranks))

    return ranking_dict


def merge_rank_dicts(rank_dicts):
    overall_ranking = dict(zip(rank_dicts[0].keys(), [[] for _ in range(len(rank_dicts[0]))]))
    for idict in range(len(rank_dicts)):
        for key in rank_dicts[idict].keys():
            overall_ranking[key].append(rank_dicts[idict][key])

    for key in overall_ranking.keys():
        overall_ranking[key] = np.mean(overall_ranking[key])

    return rank_genes(overall_ranking)#, invert_sorting=False)


def compare_highest_ranking_genes(re_path="/homes/gws/sdorkenw/rrna/data/rsubread_aligned_v2/SRR89124*fpkm",
                                  n=25):
    paths = glob.glob(re_path)

    name_list = []
    file_dict = {}
    rank_dicts = []
    for path in paths:
        name = os.path.basename(path).split("_")[0]
        name_list.append(name)
        file_dict[name] = read_values_from_file(path)
        rank_dicts.append(rank_genes(file_dict[name]))

    overall_ranking = merge_rank_dicts(rank_dicts)
    rank_to_gene = dict(zip(overall_ranking.values(), overall_ranking.keys()))

    comp = []
    gene_names = []
    for irank in range(n):
        gene_values = []
        gene_name = rank_to_gene[len(rank_to_gene)-1-irank]
        gene_names.append(gene_name)
        for path in paths:
            name = os.path.basename(path).split("_")[0]
            gene_values.append(file_dict[name][gene_name])
        comp.append(np.array(gene_values))

    comp = np.array(comp)
    return comp, name_list, gene_names


def analyse_highest_ranking_genes(n=100, plot_path=None):
    expressions, sample_names, gene_names = compare_highest_ranking_genes(n=n)
    order = np.argsort(sample_names)
    expressions = expressions[:, order]
    sample_names = np.array(sample_names)[order]

    rankings = []
    for iset in range(2):
        variations = np.std(expressions[:, 3*iset:3*(1+iset)], axis=1) / \
                     np.mean(expressions[:, 3*iset:3*(1+iset)], axis=1)

        # interesting_ids = np.argsort(variations)[:int(1*n):-1]
        interesting_ids = np.arange(n)
        interesting_genes = []
        for this_id in interesting_ids:
            interesting_genes.append(expressions[this_id, 3*iset:3*(1+iset)])

        r = np.argsort(interesting_genes, axis=1)
        print np.mean(r, axis=0), np.std(r, axis=0)
        rankings.append(r)

    if plot_path is not None:
        data = []
        for ir in range(len(rankings)):
            data += np.swapaxes(rankings[ir], 0, 1).tolist()
        data = np.array(data)
        pu.mean_std(data, sample_names, plot_path, xlabel=None,
                    ylabel="position", y_max=2)

    return rankings


def get_multiple_ribo_ratios(re_bam_path="/homes/gws/sdorkenw/rrna/data/rsubread_aligned_v2/SRR89124*_1_s.bam",
                             ribo_path="/homes/gws/sdorkenw/rrna/data/ref_genomes/rrna_hg38.gtf"):
    bam_paths = glob.glob(re_bam_path)

    ratios = {}
    for path in bam_paths:
        name = os.path.basename(path).split("_")[0]
        ratios[name] = get_ribo_ratio(path, ribo_path)

    return ratios


def get_ribo_ratio(bam_path="/homes/gws/sdorkenw/rrna/data/rsubread_aligned_v2/SRR891242_1_s.bam",
                   ribo_path="/homes/gws/sdorkenw/rrna/data/ref_genomes/rrna_hg38.gtf",
                   from_star=False):
    ribo_db = gffutils.create_db(ribo_path, dbfn=":memory:")
    ribo_features = list(ribo_db.all_features())

    pfile = pysam.AlignmentFile(bam_path, "rb")
    read_cnt = 0

    # For star aligned bam files
    cm_ids = ["CM000%d.2" % (663 + i) for i in range(24)]
    cm_mapper = {}
    for icm_id in range(1, 1 + len(cm_ids)):
        if icm_id < 23:
            cm_mapper["chr%d" % icm_id] = cm_ids[icm_id - 1]
        elif icm_id == 23:
            cm_mapper["chrX"] = cm_ids[icm_id - 1]
        else:
            cm_mapper["chrY"] = cm_ids[icm_id - 1]


    for feat in ribo_features:
        if feat.chrom in cm_mapper and feat.featuretype == "exon":
            if from_star:
                chrom = cm_mapper[feat.chrom]
            else:
                chrom = feat.chrom

            for _ in pfile.fetch(chrom, feat.start, feat.end):
                read_cnt += 1

    return read_cnt / float(pfile.mapped)


def correlate_gene_expression_depletion(n=25,
                                        re_value_path="/homes/gws/sdorkenw/rrna/data/rsubread_aligned_v2/SRR89124*fpkm",
                                        re_bam_path="/homes/gws/sdorkenw/rrna/data/rsubread_aligned_v2/SRR89124*_1_s.bam",
                                        ribo_path="/homes/gws/sdorkenw/rrna/data/ref_genomes/rrna_hg38.gtf"):
    gene_comp, sample_names, gene_names = compare_highest_ranking_genes(re_value_path, n=n)

    ribo_ratios = get_multiple_ribo_ratios(re_bam_path, ribo_path)
    ribo_ratio_list = [ribo_ratios[name] for name in sample_names]

    for igene in range(n):
        pc, pr = scipy.stats.pearsonr(ribo_ratio_list, gene_comp[igene])
        sc, sr = scipy.stats.spearmanr(ribo_ratio_list, gene_comp[igene])

        print "Rank %d: Pearson: %.4f, %.4f, Spearman: %.4f, %.4f" % (igene, pc, pr, sc, sr)


def fpkm_ribo_region_comp(re_path="/homes/gws/sdorkenw/rrna/data/rsubread_aligned_v2/SRR89124*fpkm"):
    paths = glob.glob(re_path)

    fpkm_sums = []
    for path in paths:
        fpkm_sums.append([os.path.basename(path).split("_")[0], fpkm_ribo_region(path)])

    print fpkm_sums


def fpkm_ribo_region(path):
    gene_names = ["MT-RNR1", "MT-RNR2"]

    for i in range(1, 18):
        gene_names.append("RNA5S%d" % i)

    for i in range(1, 6):
        gene_names.append("RNA18S%d" % i)

    for i in range(1, 6):
        gene_names.append("RNA28S%d" % i)

    for i in range(1, 6):
        gene_names.append("RNA-8S%d" % i)

    fpkm_dict = read_values_from_file(path)

    fpkm_sum = 0
    for name in gene_names:
        if name in  fpkm_dict:
            fpkm_sum += fpkm_dict[name]
        else:
            print "Could not find %s" % name
    return fpkm_sum