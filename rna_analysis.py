import glob
import numpy as np
import re


def read_values_from_file(path):
    gene_dict = {}

    with open(path) as f:
        for line in f.readlines():
            gene_name, value = line.split("\t")
            gene_dict[gene_name] = float(value)

    return gene_dict


def rank_genes(gene_dict):
    keys = gene_dict.keys()
    values = gene_dict.values()

    order = np.argsort(values)
    ranks = np.empty(len(order), np.int)
    ranks[order] = np.arange(len(order))
    ranking_dict = dict(zip(keys, ranks))

    return ranking_dict


def merge_rank_dicts(rank_dicts):
    overall_ranking


def compare_highest_ranking_genes(re_path):
    paths = glob.glob(re_path)

    file_dict = {}
    rank_dict = {}
    for path in paths:
        file_dict[path] = read_values_from_file(path)
        rank_dict[path] = rank_genes(file_dict[path])

