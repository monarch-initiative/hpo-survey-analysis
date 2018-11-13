import numpy as np
import argparse
from scipy.cluster.hierarchy import linkage, fcluster
from statistics import mean, median
from typing import Dict
import logging
from scipy.spatial.distance import squareform
from phenom.utils import owl_utils
from phenom.similarity.semantic_sim import SemanticSim
from rdflib import Graph
import csv

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
    parser.add_argument('--annotations', '-a', type=str, required=True,
                        help='Cached gold standard disease phenotype annotations')
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
    disease2phen = {}

    for line in ic_fh.readlines():
        hpo_id, ic = line.rstrip("\n").split("\t")
        ic_map[hpo_id] = float(ic)

    ic_fh.close()

    mondo_graph = Graph()

    #logger.info("loading mondo")
    # Previous cache made with 2018-08-03 version of mondo
    #mondo_graph.parse("/path/to/git/mondo-2018-08-03/src/ontology/reasoned.owl", format='xml')
    #root = "MONDO:0000001"

    logger.info("loading hpo")
    root = "HP:0000118"
    hpo = Graph()
    hpo.parse("http://purl.obolibrary.org/obo/hp.owl", format='xml')
    sem_sim = SemanticSim(hpo, root, ic_map)

    with open(args.annotations, 'r') as cache_file:
        reader = csv.reader(cache_file, delimiter='\t', quotechar='\"')
        for row in reader:
            if row[0].startswith('#'): continue
            (mondo_id, phenotype_id) = row[0:2]
            if mondo_id in disease2phen:
                disease2phen[mondo_id].append(phenotype_id)
            else:
                disease2phen[mondo_id] = [phenotype_id]

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

    for dist in np.linspace(50, 180, 400):
        clusters = fcluster(Z, dist, 'distance')

    #for meth in clust_meth.keys():
    #    Z = linkage(squareform(matrix), meth)
    #    clusters = fcluster(Z, 500, 'maxclust')

        cluster_map = {}
        simgic_list = []
        gic_times_cluster = []
        cluster_sizes = []
        sem_jac_list = []

        for disease_id, cluster_id in zip(labels, clusters):
            try:
                cluster_map[cluster_id].append(disease_id)
            except KeyError:
                cluster_map[cluster_id] = [disease_id]

        for cluster_id, diseases in cluster_map.items():
            cluster_sizes.append(len(diseases))
            """
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
            #mica = max([ic_map[d] for d in common_ancestors])
            #mica_list.append(mica)
            """
            profile_list = []
            for disease in diseases:
                profile_list.append(disease2phen[disease])
            sim_gic = sem_sim.groupwise_sim_gic(profile_list)
            simgic_list.append(sim_gic)
            gic_times_cluster.append(sim_gic * len(diseases))
            sim_jaccard = sem_sim.groupwise_jaccard(profile_list)
            sem_jac_list.append(sim_jaccard)


        singleton_count = sum([len(v) for k, v in cluster_map.items() if len(v) == 1])
        sizes = [len(v) for k, v in cluster_map.items()]

        #output.write("{:.4f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
        output.write("{:.4f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            #clust_meth[meth],
            dist,
            'ward',
            mean(simgic_list),
            median(simgic_list),
            len(cluster_map.keys()),
            mean(sizes),
            median(sizes),
            singleton_count,
            mean(sem_jac_list),
            median(sem_jac_list)
        ))


if __name__ == "__main__":
    main()
