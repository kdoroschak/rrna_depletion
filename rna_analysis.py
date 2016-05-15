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


def compare_highest_ranking_genes(re_path, n=25):
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


def get_multiple_ribo_ratios(re_bam_path="/homes/gws/sdorkenw/rrna/data/rsubread_aligned_v2/SRR89124*_1_s.bam",
                             ribo_path="/homes/gws/sdorkenw/rrna/data/ref_genomes/rrna_hg38.gtf"):
    bam_paths = glob.glob(re_bam_path)

    ratios = {}
    for path in bam_paths:
        name = os.path.basename(path).split("_")[0]
        ratios[name] = get_ribo_ratio(path, ribo_path)

    return ratios


def get_ribo_ratio(bam_path="/homes/gws/sdorkenw/rrna/data/rsubread_aligned_v2/SRR891242_1_s.bam",
                   ribo_path="/homes/gws/sdorkenw/rrna/data/ref_genomes/rrna_hg38.gtf"):
    ribo_db = gffutils.create_db(ribo_path, dbfn=":memory:")
    ribo_features = list(ribo_db.all_features())

    pfile = pysam.AlignmentFile(bam_path, "rb")
    read_cnt = 0

    for feat in ribo_features:
        if not "_" in feat.chrom and feat.featuretype == "exon":
            for _ in pfile.fetch(feat.chrom, feat.start, feat.end):
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

