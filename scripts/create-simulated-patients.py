from typing import Dict, Set, FrozenSet
from phenom.utils.owl_utils import get_closure
import argparse
import logging
from rdflib import Graph, RDFS
import random

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
hpo.parse("../data/hp.owl", format='xml')

top_phenotypes = {
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


def simulate_from_derived(
        pheno_profile: Set[str],
        pheno_subset: Set[str],
        graph: Graph,
        root: str,
        ic_values: Dict[str, float],
        filter_out: Set[str]) -> FrozenSet[str]:
    """
    Add imprecision and noise to profile
    20% omit phenotype - omissions
    10% use closest parent - imprecision
    30% add random phenotype - noise, min 1
    :return: FrozenSet[str] - set of phenotype curies
    """
    omission_rate = .2
    imprecision_rate = .1
    noise_rate = .3

    phenotypes = list(pheno_profile)
    profile_size = len(phenotypes)

    # Remove 20 percent of phenotypes
    count_to_remove = round(profile_size * omission_rate)
    indices = (random.sample(range(0, profile_size - count_to_remove), count_to_remove))
    for i in indices:
        del phenotypes[i]

    # mutate 10 percent to closest parent
    count_to_mutate = round(profile_size * imprecision_rate)
    #indices = (random.sample(range(0, len(list(mutatable))), count_to_mutate))
    random.shuffle(phenotypes)
    counter = 0
    for idx, pheno in enumerate(phenotypes):
        if counter == count_to_mutate:
            break
        parents = get_closure(graph, pheno, RDFS['subClassOf'], root, False)
        lay_overlap = parents.intersection(pheno_subset)
        if len(lay_overlap) == 0:
            continue
        max_ic = max([ic_values[parent] for parent in lay_overlap])
        mica = ''
        for pheno in lay_overlap:
            if ic_values[pheno] == max_ic:
                mica = pheno

        phenotypes[idx] = mica
        counter += 1

    if counter != count_to_mutate:
        logger.info("Could not mutate profile derived from {}".format(disease))

    # add random phenotype(s)
    # Filter out phenotypes from filter_out set
    phenos_to_select = pheno_subset.difference(filter_out) \
                                   .difference(phenotypes) \
                                   .difference(pheno_profile)
    if len(list(phenos_to_select)) == 0:
        logger.warn("No phenptypes to select for "
                    "profile derived from {}".format(disease))
    comissions = round(profile_size * noise_rate)
    noise_count = 1 if comissions == 0 else comissions
    for i in range(noise_count):
        random_pheno = random.choice(list(phenos_to_select))
        phenotypes.append(random_pheno)
        phenos_to_select = phenos_to_select - set(random_pheno)

    return frozenset(phenotypes)


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
                filter_out = top_phenotypes
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
