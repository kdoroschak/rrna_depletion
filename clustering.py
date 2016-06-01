import cPickle as pkl
import glob
import h5py
import numpy as np
import os
import re
from scipy.spatial import distance
from sklearn.cluster import KMeans
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA
import sklearn.preprocessing

def iterate_over_n_cluster(data, n_clusters):
    bics = []
    data = sklearn.preprocessing.scale(data)
    for n_cluster in n_clusters:
        kmeans = calc_kmeans(data, n_cluster=n_cluster, init="k-means++", save=False)
        bics.append(compute_bic(kmeans[1], data))
        print n_cluster, bics[-1]
    return bics


def calc_kmeans(data, n_cluster=4, init="k-means++", save=True, path_suffix=""):
    # data = np.swapaxes(data, 0, 1)

    if init == "pca":
        pca = PCA(n_components=n_cluster).fit(data)
        est = KMeans(init=pca.components_, n_clusters=n_cluster, n_init=1, n_jobs=20)
    elif init == "k-means++":
        est = KMeans(init='k-means++', n_clusters=n_cluster, n_init=4, n_jobs=20)
    elif init == "random":
       est = KMeans(init="random", n_clusters=n_cluster, n_init=16, n_jobs=20)
    else:
        raise Exception("init not known")

    fit = est.fit(data)
    prediction = est.predict(data)

    if save:
        with open("/homes/gws/sdorkenw/rrna/data/clusterings/fit" + path_suffix + ".pkl", "w") as f:
            pkl.dump(fit, f)
        with open("/homes/gws/sdorkenw/rrna/data/clusterings/prediction" + path_suffix + ".pkl", "w") as f:
            pkl.dump(prediction, f)
    else:
        return est, fit, prediction


def compute_bic(kmeans, X):
    """ http://stats.stackexchange.com/questions/90769/using-bic-to-estimate-the-number-of-k-in-kmeans
    Computes the BIC metric for a given clusters

    Parameters:
    -----------------------------------------
    kmeans:  List of clustering object from scikit learn

    X     :  multidimension np array of data points

    Returns:
    -----------------------------------------
    BIC value
    """
    # assign centers and labels
    centers = [kmeans.cluster_centers_]
    labels = kmeans.labels_
    # number of clusters
    m = kmeans.n_clusters
    # size of the clusters
    n = np.bincount(labels)
    # size of data set
    N, d = X.shape

    # compute variance for all clusters beforehand
    cl_var = (1.0 / (N - m) / d) * sum([sum(distance.cdist(X[np.where(labels == i)],
                                                           [centers[0][i]],
                                                           'euclidean')**2)
                                        for i in range(m)])[0]

    const_term = 0.5 * m * np.log(N) * (d+1)

    bic = np.sum([n[i] * np.log(n[i]) -
                  n[i] * np.log(N) -
                  ((n[i] * d) / 2.) * np.log(2*np.pi*cl_var) -
                  ((n[i] - 1) * d / 2.) for i in range(m)]) - const_term

    return bic