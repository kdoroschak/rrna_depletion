from sklearn.neighbors import NearestNeighbors
import cPickle as pkl


class RiboNN(object):
    def __init__(self, data):
        self.nn = None
        self.data = data

    def init_with_data(self, data):
        neigh = NearestNeighbors(5, 5, metric="hamming", n_jobs=20)
        self.nn = neigh.fit(data)
