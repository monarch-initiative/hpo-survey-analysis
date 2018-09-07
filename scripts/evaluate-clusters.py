import numpy as np
import argparse
from scipy.cluster.hierarchy import linkage, fcluster
from statistics import mean, median
from typing import Dict
import logging
from scipy.spatial.distance import squareform
from phenom.utils import owl_utils
from rdflib import Graph

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    """
    Cluster distance matrix with scipy.cluster.hierarchy
    and evaluate clusters with information content of
    MONDO classes
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

    ic_fh = open(args.ic_cache, 'r')
    output = open(args.output, 'w')
    output.write("#distance\tlinkage\tmean_mica\tmedian_mica\t"
                 "num_clusters\tmean_mem\tmedian_mem\tsingletons\n")

    ic_map: Dict[str, float] = {}

    for line in ic_fh.readlines():
        hpo_id, ic = line.rstrip("\n").split("\t")
        ic_map[hpo_id] = float(ic)

    ic_fh.close()

    mondo_graph = Graph()

    logger.info("loading mondo")
    # Previous cache made with 2018-08-03 version of mondo
    mondo_graph.parse("/home/kshefchek/git/mondo-2018-08-03/src/ontology/reasoned.owl", format='xml')
    root = "MONDO:0000001"

    mondo_skip = {
        'MONDO:0023807',
        'MONDO:0000559',
        'MONDO:0009117',
        'MONDO:0016961',
        'MONDO:0011750',
        'MONDO:0017180',
    }
    clust_meth = {
        'ward': 'ward',
        'average': 'UPGMA',
        'weighted': 'WPGMA',
        'centroid': 'centroid',
        'complete': 'complete',
        'median': 'median'
    }

    logger.info("clustering")
    Z = linkage(squareform(matrix), 'ward')

    for dist in np.linspace(.35, .5, 91):
        clusters = fcluster(Z, dist, 'distance')

    #for meth in clust_meth.keys():
    #    Z = linkage(squareform(matrix), meth)
    #    clusters = fcluster(Z, 500, 'maxclust')

        cluster_map = {}
        mica_list = []
        cluster_sizes = []
        unclustered = 0
        poor = 0
        medium = 0
        good = 0
        excellent = 0

        for disease_id, cluster_id in zip(labels, clusters):
            try:
                cluster_map[cluster_id].append(disease_id)
            except KeyError:
                cluster_map[cluster_id] = [disease_id]

        for cluster_id, diseases in cluster_map.items():
            cluster_sizes.append(len(diseases))
            if len(set(diseases).intersection(mondo_skip)) > 0:
                if len(diseases) == len(set(diseases).intersection(mondo_skip)):
                    logger.warning("Cannot evaluate cluster")
                    diseases = set(diseases)
                else:
                    diseases = set(diseases) - mondo_skip
            common_ancestors = set()
            is_first = True
            for disease in diseases:
                if is_first:
                    common_ancestors = owl_utils.get_closure(mondo_graph, disease, root=root)
                    is_first = False
                else:
                    common_ancestors = common_ancestors.intersection(
                        owl_utils.get_closure(mondo_graph, disease, root=root)
                    )
            mica = max([ic_map[d] for d in common_ancestors])
            mica_list.append(mica)
            if mica <= .001:
                unclustered += len(cluster_map[cluster_id])
            elif .001 < mica <= 3:
                poor += len(cluster_map[cluster_id])
            elif 3 < mica <= 6:
                medium += len(cluster_map[cluster_id])
            elif 6 < mica <= 9:
                good += len(cluster_map[cluster_id])
            else:
                excellent += len(cluster_map[cluster_id])

        singleton_count = sum([len(v) for k, v in cluster_map.items() if len(v) == 1])
        sizes = [len(v) for k, v in cluster_map.items()]

        #output.write("{:.4f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
        output.write("{:.4f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            #clust_meth[meth],
            dist,
            'ward',
            mean(mica_list),
            median(mica_list),
            len(cluster_map.keys()),
            mean(sizes),
            median(sizes),
            singleton_count,
            unclustered,
            poor,
            medium,
            good,
            excellent
        ))


if __name__ == "__main__":
    main()
