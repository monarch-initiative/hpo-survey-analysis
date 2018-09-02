from phenom.similarity.semantic_sim import SemanticSim
from phenom.similarity.semantic_dist import SemanticDist
import argparse
import logging
import csv
from rdflib import Graph
from typing import Dict, List, Tuple, Any
from itertools import combinations
import multiprocessing
from multiprocessing import Process, Queue

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description='Given a subset of HPO terms and diseases, generates '
                    'derived annotations from the HPO disease to phenotype '
                    'annotations')
    parser.add_argument('--phenotypes', '-ph', type=str, required=False)
    parser.add_argument('--diseases', '-d', type=str, required=True)
    parser.add_argument('--ic_cache', '-ic', type=str, required=True)
    parser.add_argument('--annotations', '-a', type=str, required=True,
                    help='Cached gold standard disease phenotype annotations')
    parser.add_argument('--processes', '-p', type=int, required=False,
                    default=int(multiprocessing.cpu_count()/2),
                    help='Number of processes to spawn')
    parser.add_argument('--output', '-o', type=str, required=False,
                        help='Location of output file', default="./matrix.csv")

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

    dist_matrix = [[0 for k in range(len(diseases))] for i in range(len(diseases))]
    ic_map: Dict[str, float] = {}
    disease2phen: Dict[str, List[str]] = {}
    disease_map: Dict[str, int] = {}

    # create dictionary of diseases and their index
    index = 0
    for disease in diseases:
        disease_map[disease] = index
        index += 1

    # Create list of combinations of diseases
    disease_combos = list(combinations(diseases, 2))

    for line in ic_fh.readlines():
        hpo_id, ic = line.rstrip("\n").split("\t")
        ic_map[hpo_id] = float(ic)

    with open(args.annotations, 'r') as cache_file:
        reader = csv.reader(cache_file, delimiter='\t', quotechar='\"')
        for row in reader:
            if row[0].startswith('#'): continue
            (mondo_id, phenotype_id) = row[0:2]
            if mondo_id in disease2phen:
                disease2phen[mondo_id].append(phenotype_id)
            else:
                disease2phen[mondo_id] = [phenotype_id]

    procs = []
    queue = Queue()  # create a queue object
    result_list = []

    # Split into chunks depending on args.processes
    for chunk in [disease_combos[i::args.processes] for i in range(args.processes)]:
        proc = Process(target=get_resnik_sim,
                       args=(chunk, disease2phen,hpo, root, ic_map, disease_map, queue))
        proc.start()
        procs.append(proc)

    for i in range(args.processes):
        result_list.extend(queue.get())

    for proc in procs:
        proc.join()

    for dis_a_index, dis_b_index, score in result_list:
        dist_matrix[dis_a_index][dis_b_index] = score
        dist_matrix[dis_b_index][dis_a_index] = score

    for row in dist_matrix:
        csv_writer.writerow(row)


def get_resnik_sim(combos, disease2phen, graph, root, ic_map, coordinates, queue):
    #sem_sim = SemanticSim(graph, root, ic_map)
    sem_dist = SemanticDist(graph, root, ic_map)
    result_list: List[Tuple[int,int,Any]] = []
    for index, combo in enumerate(combos):
        disease_a, disease_b = combo
        disease_a_profile = disease2phen[disease_a]
        disease_b_profile = disease2phen[disease_b]
        #score = 1 - sem_sim.resnik_sim(
        #    disease_a_profile, disease_b_profile, is_normalized=True, is_symmetric=True)
        #score = 1 - sem_sim.cosine_sim(
        #    disease_a_profile, disease_b_profile, ic_weighted=True)
        score = sem_dist.euclidean_matrix(disease_a_profile, disease_b_profile)
        if score == 1 or score == 0:
            result_list.append((coordinates[disease_a],
                               coordinates[disease_b],
                               int(score))
            )
        else:
            result_list.append((coordinates[disease_a],
                               coordinates[disease_b],
                               "{:.4f}".format(score))
            )

        if index % 100000 == 0:
            logger.info("Processed {} combinations out of {}".format(index, len(combos)))

    queue.put(result_list)


if __name__ == "__main__":
    main()
