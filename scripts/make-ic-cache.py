from typing import Dict, Set
import csv
import requests
from contextlib import closing
import codecs
from phenom import monarch
import logging
import argparse
from argparse import ArgumentError

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# globals and constants
HPO_ONTOLOGY = "http://purl.obolibrary.org/obo/hp.owl"

# contains manual and semi-automated annotations created by the HPO-team.
# These are annotations of OMIM-, Orphanet-, and DECIPHER-entries
HPO_DATA = "http://compbio.charite.de/jenkins/job/hpo.annotations/" \
           "lastStableBuild/artifact/misc/phenotype_annotation.tab"


def main():

    parser = argparse.ArgumentParser(
        description='Generate information content for each HPO class using the '
                    'HPO phenotype annotation file ')
    parser.add_argument('--mondo_cache', '-m', type=str, required=False,
                        help='Cached 2 column disease phenotype tsv')
    parser.add_argument('--mondo_output', '-mo', type=str, required=False,
                        help='path to mondo 2 column file (if no cache)')
    parser.add_argument('--output', '-o', type=str, required=False,
                        help='Location of output file', default="./output.tsv")
    args = parser.parse_args()

    if not args.mondo_cache and not args.mondo_output:
        raise ArgumentError(
            args.mondo_cache, 'One of --mondo_cache or '
                              '--mondo_output must be provided')

    # i/o
    output_file = open(args.output, 'w')

    d2p_map = dict()
    if not args.mondo_cache:
        mondo_cache = open(args.mondo_output, 'w')
        d2p_map = _process_hpo_data(HPO_DATA)
        for disease, pheno_set in d2p_map.items():
            for pheno in pheno_set:
                mondo_cache.write("{}\t{}\n".format(disease, pheno))
    else:
        with open(args.mondo_cache, 'r') as cache_file:
            pass


def _process_hpo_data(url: str) -> Dict[str, Set[str]]:

    disease2pheno = dict()

    # https://stackoverflow.com/a/38677619
    with closing(requests.get(url, stream=True)) as r:
        counter = 0
        reader = csv.reader(codecs.iterdecode(r.iter_lines(), 'utf-8'), delimiter='\t', quotechar='\"')
        for row in reader:
            (db, num, name, qual, pheno_id) = row[0:5]
            disease_id = "{}:{}".format(db, num)

            mondo_node = monarch.get_clique_leader(disease_id)
            mondo_id = mondo_node['id']

            if mondo_id in disease2pheno:
                disease2pheno[mondo_id].add(pheno_id)
            else:
                disease2pheno[mondo_id] = {pheno_id}

            counter += 1
            if counter % 10000 == 0:
                print("processed {} rows".format(counter))

    return disease2pheno


if __name__=="__main__":
    main()