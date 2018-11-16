from phenom import monarch
from phenom.utils import owl_utils
import argparse
import logging
import csv
from rdflib import Graph, RDFS
from typing import Dict, List

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description='Given a subset of HPO terms and diseases, generates '
                'derived annotations from the HPO disease to phenotype annotations')
parser.add_argument('--phenotypes', '-p', type=str, required=True)
parser.add_argument('--diseases', '-d', type=str, required=True)
parser.add_argument('--ic_cache', '-ic', type=str, required=True)
parser.add_argument('--annotations', '-a', type=str, required=False,
                    help='Cached gold standard disease phenotype annotations')
parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./derived-cache.tsv")

args = parser.parse_args()

root = "HP:0000118"
hpo = Graph()
hpo.parse("../data/hp.owl", format='xml')

# Dictionaries
ic_map: Dict[str, float] = {}
gold_standard: Dict[str, List[str]] = {}

# I/O
pheno_fh = open(args.phenotypes, 'r')
disease_fh = open(args.diseases, 'r')
ic_fh = open(args.ic_cache, 'r')
output = open(args.output, 'w')
lay_terms = set(pheno_fh.read().splitlines())
diseases = disease_fh.read().splitlines()


for line in ic_fh.readlines():
    hpo_id, ic = line.rstrip("\n").split("\t")
    ic_map[hpo_id] = float(ic)

if args.annotations:
    with open(args.annotations, 'r') as cache_file:
        reader = csv.reader(cache_file, delimiter='\t', quotechar='\"')
        for row in reader:
            if row[0].startswith('#'): continue
            (mondo_id, phenotype_id) = row[0:2]
            if mondo_id in gold_standard:
                gold_standard[mondo_id].append(phenotype_id)
            else:
                gold_standard[mondo_id] = [phenotype_id]

# Load from solr (note these may be from an older version of HPOA)
else:
    for mondo in diseases:
        # Get phenotypes
        pheno_profile, mondo_label = monarch.get_direct_phenotypes(mondo)
        gold_standard[mondo] = pheno_profile


for mondo in diseases:
    # If list is not mondo
    # clique_leader = monarch.get_clique_leader(disease)
    # mondo = clique_leader['id']
    # mondo_label = clique_leader['label']

    # Get phenotypes
    gold_profile = set(gold_standard[mondo])
    derived_profile = gold_profile.intersection(lay_terms)
    non_lay_terms = gold_profile - derived_profile

    for phenotype in non_lay_terms:
        parents = owl_utils.get_closure(hpo, phenotype, RDFS['subClassOf'], root)
        lay_overlap = parents.intersection(lay_terms)
        if len(lay_overlap) == 0:
            continue
        max_ic = max([ic_map[parent] for parent in lay_overlap])
        mica = ''
        for pheno in lay_overlap:
            if ic_map[pheno] == max_ic:
                mica = pheno

        derived_profile.add(mica)

    for phenotype in derived_profile:
        output.write("{}\t{}\n".format(mondo, phenotype))
