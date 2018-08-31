from phenom import monarch
from phenom.similarity.semantic_sim import SemanticSim
from phenom.utils import owl_utils
import argparse
import logging
import csv
from rdflib import Graph
from typing import Dict, List
from collections import deque
import multiprocessing
from multiprocessing import Process, Manager

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_score(query_profile, diseases, scores):
    sem_sim = SemanticSim(hpo, root, ic_map)
    for disease in diseases:
        pheno_profile = gold_standard[disease]
        score = 1 - sem_sim.resnik_sim(
            query_profile, pheno_profile, is_normalized=True, is_symmetric=True)
        scores.append((disease, score))
    return scores


parser = argparse.ArgumentParser(
    description='Given a subset of HPO terms and diseases, generates '
                'derived annotations from the HPO disease to phenotype annotations')
parser.add_argument('--phenotypes', '-ph', type=str, required=False)
parser.add_argument('--diseases', '-d', type=str, required=True)
parser.add_argument('--ic_cache', '-ic', type=str, required=True)
parser.add_argument('--annotations', '-a', type=str, required=True,
                    help='Cached gold standard disease phenotype annotations')
parser.add_argument('--processes', '-p', type=int, required=False,
                    default=int(multiprocessing.cpu_count()/2),
help='Number of processes to spawn')

parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./matrix.tsv")

args = parser.parse_args()


root = "HP:0000118"
hpo = Graph()
hpo.parse("http://purl.obolibrary.org/obo/hp.owl", format='xml')

# I/O
disease_fh = open(args.diseases, 'r')
ic_fh = open(args.ic_cache, 'r')
output = open(args.output, 'w')
csv_writer = csv.writer(output, delimiter=',')
diseases = disease_fh.read().splitlines()

ic_map: Dict[str, float] = {}
gold_standard: Dict[str, List[str]] = {}

for line in ic_fh.readlines():
    hpo_id, ic = line.rstrip("\n").split("\t")
    ic_map[hpo_id] = float(ic)

if args.phenotypes:
    pheno_fh = open(args.phenotypes, 'r')
    phenotype_terms = set(pheno_fh.read().splitlines())
else:
    phenotype_terms = owl_utils.get_descendants(hpo, root)

# memoize closures
# doesn't quite work with multiprocessing
# for terms in phenotype_terms:
#    owl_utils.get_closure(hpo, terms, root=root)

with open(args.annotations, 'r') as cache_file:
    reader = csv.reader(cache_file, delimiter='\t', quotechar='\"')
    for row in reader:
        if row[0].startswith('#'): continue
        (mondo_id, phenotype_id) = row[0:2]
        if mondo_id in gold_standard:
            gold_standard[mondo_id].append(phenotype_id)
        else:
            gold_standard[mondo_id] = [phenotype_id]

gs = gold_standard['MONDO:0000044']

is_first_run = True

for disease_index in range(len(diseases)):

    manager = Manager()
    scores = manager.list()
    procs = []
    disease_score: Dict[str, float] = {}

    if not is_first_run:
        del diseases[0]
    else:
        is_first_run = False

    matrix_row = [0 for i in range(disease_index)]
    query_profile = gold_standard[diseases[disease_index]]

    # Split into chunks depending on args.processes
    for chunk in [diseases[i::args.processes] for i in range(args.processes)]:
        proc = Process(target=get_score, args=(query_profile, chunk, scores))
        proc.start()
        procs.append(proc)

    for proc in procs:
        proc.join()

    for disease, score in scores:
        disease_score[disease] = score

    for disease in diseases:
        matrix_row.append(disease_score[disease])

    csv_writer.writerow(matrix_row)
    if disease_index % 100 == 0:
        logger.info("Processed {} rows".format(disease_index))

"""
scores = []
for disease in diseases:
    pheno_profile = gold_standard[disease]
    score = sem_sim.resnik_sim(gs, pheno_profile, is_normalized=True, is_symmetric=True)
    scores.append(score)
"""

