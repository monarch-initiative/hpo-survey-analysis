import numpy as np
import argparse
from scipy.cluster.hierarchy import linkage, fcluster
from statistics import mean, median
from typing import Dict
import logging
from scipy.spatial.distance import squareform
from phenom.utils import owl_utils
from rdflib import Graph
from phenom import monarch


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    """
    Cluster and iterate over each cluster to find the best
    disease group that subsumes the cluster
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

    cluster_map = {}

    logger.info("clustering")
    Z = linkage(squareform(matrix), 'ward')
    # cosine weighted = 2631
    # resnik = 2453
    # euclidean = 525
    clusters = fcluster(Z, 525, 'maxclust')

    for disease_id, cluster_id in zip(labels, clusters):
        try:
            cluster_map[cluster_id].append(disease_id)
        except KeyError:
            cluster_map[cluster_id] = [disease_id]

    for cluster_id, diseases in cluster_map.items():
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
        for dis in common_ancestors:
            if ic_map[dis] == mica:
                mica_id = dis
        label = monarch.get_label(mica_id)
        # Number of subclasses for dis
        subclass_count = len(owl_utils.get_closure(mondo_graph, mica_id, negative=True, reflexive=False))
        output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            cluster_id,
            len(cluster_map[cluster_id]),
            mica_id,
            label,
            mica,
            subclass_count,
            "|".join(cluster_map[cluster_id])
        ))


if __name__ == "__main__":
    main()