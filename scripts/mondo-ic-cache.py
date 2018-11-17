from typing import Dict
import csv
from phenom.utils import owl_utils
from phenom.math import math_utils
from rdflib import Graph
import logging
import argparse


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description='Generate information content for each HPO class using the '
                    'HPO phenotype annotation file ')
    parser.add_argument('--mondo_cache', '-m', type=str, required=True,
                        help='Cached 2 column disease phenotype tsv')
    parser.add_argument('--output', '-o', type=str, required=False,
                        help='Location of output file', default="./mondo-ic-cache.tsv")
    args = parser.parse_args()

    # i/o
    output_file = open(args.output, 'w')

    explicit_annotations = 1
    disease_annotations: Dict[str, int] = {}

    mondo_graph = Graph()

    # Previous cache made with 2018-08-03 version of mondo
    logger.info("Loading MONDO")
    mondo_graph.parse("/path/to/git/mondo-2018-08-03/src/ontology/reasoned.owl", format='xml')
    root = "MONDO:0000001"

    logger.info("Getting classes")
    all_diseases = owl_utils.get_descendants(mondo_graph, root)
    disease_annotations = {disease: 0 for disease in all_diseases}

    logger.info("Seeding leaf nodes")
    # Seed leaf nodes with 1 annotation
    for leaf in owl_utils.get_leaf_nodes(mondo_graph, root):
        explicit_annotations += 1
        for disease in owl_utils.get_closure(mondo_graph, leaf, root=root):
            try:
                disease_annotations[disease] += 1
            except KeyError:
                print(disease)
                disease_annotations[disease] = 1

    logger.info("Fetching annotations")
    with open(args.mondo_cache, 'r') as cache_file:
        reader = csv.reader(cache_file, delimiter='\t', quotechar='\"')
        for row in reader:
            if row[0].startswith('#'): continue
            if not row[0].startswith('MONDO'): continue
            (mondo_id, phenotype_id) = row[0:2]
            explicit_annotations += 1
            for disease in owl_utils.get_closure(mondo_graph, mondo_id, root=root):
                try:
                    disease_annotations[disease] += 1
                except KeyError:
                    print(disease)
                    disease_annotations[disease] = 1

    logger.info("Computing IC")
    for disease, annot_count in disease_annotations.items():
        output_file.write("{}\t{}\n".format(
            disease,
            math_utils.information_content(annot_count / explicit_annotations)
        ))


if __name__ == "__main__":
    main()