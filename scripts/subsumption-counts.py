from phenom.utils.owl_utils import get_closure
from rdflib import Graph
from typing import Dict, Set
import argparse

parser = argparse.ArgumentParser(
    description='Determines by uniqueness counting the number'
                ' of other diseases that subsume it')
parser.add_argument('--derived_annotations', '-da', type=str, required=True)
parser.add_argument('--annotations', '-a', type=str, required=True,
                    help='Cached gold standard disease phenotype annotations')
parser.add_argument('--labels', '-l', type=str, required=True,
                    help='Mondo labels')
parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./subsumption-counts.tsv")
args = parser.parse_args()

# i/o
output = open(args.output, "w")
output_fields = [
    'derived from disease id',
    'disease label',
    'derived subsumed by',
    'gold subsumed by'
]
output.write("{}\n".format("\t".join(output_fields)))

hpo = Graph()
hpo.parse("../data/hp.owl", format='xml')

gold_standard: Dict[str, Set[str]] = {}
derived_profiles: Dict[str, Set[str]] = {}
mondo_diseases: Dict[str, str] = {}

if args.labels:
    with open(args.labels, 'r') as mondo_labels:
        for line in mondo_labels:
            disease_id, disease_label = line.rstrip("\n").split("\t")
            mondo_diseases[disease_id] = disease_label


def load_map_from_file(file_path: str) -> Dict[str, Set[str]]:
    profile_map: Dict[str, Set[str]] = {}
    with open(file_path, 'r') as annotations:
        for line in annotations:
            if line.startswith('#') or not line.startswith('MONDO'):
                continue
            disease, phenotype = line.rstrip("\n").split("\t")[0:2]
            try:
                profile_map[disease].add(phenotype)
            except KeyError:
                profile_map[disease] = {phenotype}
            for pheno in get_closure(hpo, phenotype, root='HP:0000118', reflexive=False):
                profile_map[disease].add(pheno)
    return profile_map


gold_standard = load_map_from_file(args.annotations)
derived_profiles = load_map_from_file(args.derived_annotations)

for disease, profile in derived_profiles.items():
    sub_count = 0
    gold_count = 0
    gold_profile = gold_standard[disease]
    for gold_disease, hpo_profile in gold_standard.items():
        if len(profile & hpo_profile) == len(profile):
            sub_count += 1
        if len(gold_profile & hpo_profile) == len(gold_profile):
            gold_count += 1
    try:
        label = mondo_diseases[disease]
    except KeyError:
        print(disease)
        label = 'obsoleted class'
    output.write("{}\t{}\t{}\t{}\n".format(disease, label, sub_count, gold_count))
