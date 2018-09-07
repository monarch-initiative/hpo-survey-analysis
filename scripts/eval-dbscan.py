import numpy as np
import argparse
from statistics import mean, median
import logging
from sklearn.cluster import DBSCAN

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    """
    Cluster distance matrix with dbscan
    """
    parser = argparse.ArgumentParser(description='description')
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file'
                             ' that contains the sim matrix as json')
    parser.add_argument('--label', '-l', type=str, required=True,
                        help='Location of id-label mapping file')
    parser.add_argument('--ic_cache', '-ic', type=str, required=True)
    parser.add_argument('--output', '-o', required=False, help='output file')
    args = parser.parse_args()

    logger.info("loading matrix")
    matrix = np.loadtxt(args.input, delimiter=",")
    labels = [line.rstrip('\n').split('\t')[0] for line in open(args.label, 'r')]

    cluster_map = {}

    db = DBSCAN(eps=.32, metric="precomputed").fit(matrix)
    singleton = -1
    for disease_id, cluster_id in zip(labels, db.labels_):

        try:
            if cluster_id == -1:
                cluster_map[singleton] = [disease_id]
                singleton -= 1
            else:
                cluster_map[cluster_id].append(disease_id)
        except KeyError:
            cluster_map[cluster_id] = [disease_id]

    singleton_count = sum([len(v) for k, v in cluster_map.items() if len(v) == 1])
    sizes = [len(v) for k, v in cluster_map.items()]

    logger.info("{} clusters".format(len(cluster_map.keys())))
    logger.info("{} singletons".format(singleton_count))
    logger.info("Avg cluster size: {}".format(mean(sizes)))
    logger.info("median cluster size: {}".format(median(sizes)))


if __name__ == "__main__":
    main()
