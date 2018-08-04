from phenom import monarch
from phenom.similarity.phenodigm import Phenodigm
from phenom.utils import owl_utils
import argparse
import logging
from rdflib import Graph, RDFS
import random

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description='Given a subset of HPO terms, generates derived annotations '
                'from the HPO disease to phenotype annotations and computes '
                'similarity and sufficiency scores for each disease')
parser.add_argument('--phenotypes', '-p', type=str, required=True)
parser.add_argument('--diseases', '-d', type=str, required=True)
parser.add_argument('--ic_cache', '-ic', type=str, required=True)
parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./output.tsv")

args = parser.parse_args()

root = "HP:0000118"
hpo = Graph()
hpo.parse("http://purl.obolibrary.org/obo/hp.owl", format='xml')


# I/O
pheno_fh = open(args.phenotypes, 'r')
disease_fh = open(args.diseases, 'r')
ic_fh = open(args.ic_cache, 'r')
ic_map = {}

for line in ic_fh.readlines():
    hpo_id, ic = line.rstrip("\n").split("\t")
    ic_map[hpo_id] = float(ic)

output = open(args.output, 'w')

pheno_list = set(pheno_fh.read().splitlines())
disease_list = disease_fh.read().splitlines()

phenodigm = Phenodigm(hpo, root, ic_map)

for mondo in disease_list:
    # Get clique leader if list is not mondo
    # clique_leader = monarch.get_clique_leader(disease)
    # mondo = clique_leader['id']
    # mondo_label = clique_leader['label']

    # Get phenotypes
    pheno_profile, mondo_label = monarch.get_direct_phenotypes(mondo)

    lay_profile = pheno_profile.intersection(pheno_list)

    disease_only = pheno_profile - lay_profile

    for phenotype in disease_only:
        parents = owl_utils.get_closure(hpo, phenotype, RDFS['subClassOf'], root)
        lay_overlap = parents.intersection(pheno_list)
        if len(lay_overlap) == 0:
            continue
        max_ic = max([ic_map[parent] for parent in lay_overlap])
        mica = ''
        for pheno in lay_overlap:
            if ic_map[pheno] == max_ic:
                mica = pheno
        lay_profile.add(mica)

    # Get annot sufficiency of whole profile
    scores = monarch.get_annotation_sufficiency_score(pheno_profile)
    disease_score = scores['scaled_score']

    # Add noise to profile
    """
    Remove 5% of annotations at random
    If profile is <= 20: 
       randomly select 10% phenotypes and move them
       one level higher in the HPO ontology

    If the profile contains 1-4 phenotypes, move one level higher in the ontology
    """
    lay_profile = list(lay_profile)
    profile_size= len(lay_profile)
    if profile_size >= 20:
        count_to_remove = round(profile_size * .05)
        indices = (random.sample(range(0, profile_size-count_to_remove), count_to_remove))
        for i in indices:
            del lay_profile[i]

    elif profile_size < 10 and profile_size != 0:
        count_to_mutate = round(profile_size * .1)
        if count_to_mutate == 0:
            count_to_mutate = 1
        indices = (random.sample(range(0, profile_size), count_to_mutate))
        for i in indices:
            sel_phenotype = lay_profile[i]
            # get parent(s) select 1st if multiple
            objs = hpo.objects(owl_utils.hp_curie2iri(sel_phenotype), RDFS['subClassOf'])
            lay_profile[i] = owl_utils.hp_iri2curie(str(list(objs)[0]))

    if (len(lay_profile) == 0):
        output.write("{}\t{}\t{}\n".format(
            mondo,
            mondo_label,
            "0",
        ))
        continue

    try:
        sim_score = phenodigm.symmetric_phenodigm(lay_profile, pheno_profile)
        #sim_score = monarch.owlsim_compare(lay_profile, pheno_profile)
    except ValueError:
        print(lay_profile)
        print(pheno_profile)
    except IndexError as e:
        print(lay_profile)
        print(pheno_profile)
        raise e

    output.write("{}\t{}\t{:0.3f}\n".format(
        mondo,
        mondo_label,
        sim_score,
    ))

