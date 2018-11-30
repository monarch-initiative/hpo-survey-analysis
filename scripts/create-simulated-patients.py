from typing import Dict, Set
import argparse
import logging
from rdflib import Graph
from phenom.utils.simulate import simulate_from_derived

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description='Given derived annotations, simulates patients '
                'by adding imprecision, noise, and omissions')
parser.add_argument('--derived', '-d', type=str, required=True)
parser.add_argument('--phenotypes', '-p', type=str, required=True)
parser.add_argument('--ic_cache', '-ic', type=str, required=True)
parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./simulated-patients.tsv")

args = parser.parse_args()

patients_per_disease = 20

phenotype_subset: Set[str] = set()
derived_profiles: Dict[str, Set[str]] = {}
ic_map: Dict[str, float] = {}


abn_phenotype = "HP:0000118"
hpo = Graph()
hpo.parse("../data/owl/hp.owl", format='xml')

top_phenotypes = {
    "HP:0000118",
    'HP:0000924',
    'HP:0000707',
    'HP:0000152',
    'HP:0001574',
    'HP:0000478',
    'HP:0001626',
    'HP:0001939',
    'HP:0000119',
    'HP:0025031',
    'HP:0002664',
    'HP:0001871',
    'HP:0002715',
    'HP:0000818',
    'HP:0003011',
    'HP:0002086',
    'HP:0000598',
    'HP:0003549',
    'HP:0001197',
    'HP:0001507',
    'HP:0000769'
}

# I/O
pheno_fh = open(args.phenotypes, 'r')
with open(args.phenotypes, 'r') as pheno_file:
    for line in pheno_file:
        phenotype_subset.add(line.rstrip("\n"))

with open(args.derived, 'r') as annotations:
    for line in annotations:
        if line.startswith('#') or not line.startswith('MONDO'):
            continue
        disease, phenotype = line.rstrip("\n").split("\t")[0:2]
        try:
            derived_profiles[disease].add(phenotype)
        except KeyError:
            derived_profiles[disease] = {phenotype}

with open(args.ic_cache, 'r') as ic_file:
    for line in ic_file:
        hpo_id, ic = line.rstrip("\n").split("\t")
        ic_map[hpo_id] = float(ic)

output = open(args.output, 'w')

logging.info("creating simulated patients")
sim_count = 0
total_patients = len(derived_profiles.keys()) * patients_per_disease
prof_w_one_pheno = 0
for disease, derived_prof in derived_profiles.items():
    if len(list(derived_prof)) == 1:
        prof_w_one_pheno += 1

    max_iterations = 200
    iterations = 0
    if sim_count % 10000 == 0:
        logger.info("Created {} patients out of {}".format(
            sim_count,
            total_patients
        ))

    simulated_patients = set()
    while iterations != max_iterations:
        simulated_patients.add(
            simulate_from_derived(
                pheno_profile = derived_prof,
                pheno_subset = phenotype_subset,
                graph = hpo,
                root = abn_phenotype,
                ic_values = ic_map,
                filter_out = top_phenotypes,
                ref_disease = disease
            )
        )
        iterations += 1
        if len(list(simulated_patients)) == patients_per_disease:
            break

    if iterations == max_iterations:
        logger.info("Could not create {} patients from {}, created {}".format(
            patients_per_disease,
            disease,
            len(list(simulated_patients))
        ))

    patient_key = 0
    for sim in simulated_patients:
        for pheno in sim:
            output.write("{0}-{1}\t{2}\t{0}\n".format(
                disease,
                patient_key,
                pheno
            ))
        patient_key += 1
    sim_count += patients_per_disease

logger.info("{} derived profiles have 1 phenotype".format(prof_w_one_pheno))
