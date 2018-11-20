"""Given a file containing synthetic patient data, generates
confusion matrix per rank for the owlsim3 bayes classifier
"""
import argparse
import logging
from typing import Dict, Set, List
import gzip
import multiprocessing
from multiprocessing import Process, Queue
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description='Create a table of true positive/neg and false pos/neg '
                'for each rank in a binary classifier')
parser.add_argument('--simulated', '-s', type=str, required=True)
parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./confusion-matrix.tsv")
parser.add_argument('--processes', '-p', type=int, required=False,
                    default=int(multiprocessing.cpu_count()/2))

args = parser.parse_args()


def create_confusion_matrix_per_rank(
        synthetic_patients: List,
        owlsim_url: str,
        sim_profiles: Dict[str, Set[str]],
        disease_synth_map: Dict[str, str],
        ranks_to_eval: int,
        queue: Queue):

    confusion_by_rank: Dict[int, List[int]] = {}

    for rank in range(1, ranks_to_eval + 1):
        # tp, fp, fn, tn
        confusion_by_rank[rank] = [0, 0, 0, 0]

    counter = 0
    total = len(synthetic_patients)

    for synth_patient in synthetic_patients:

        if counter % 1000 == 0:
            logging.info("processed {} patients out of {}".format(counter, total))
        counter += 1

        # Useful for testing
        if counter == 500:
            break

        params = {
            'id': sim_profiles[synth_patient],
            'limit': 8000
        }
        sim_req = requests.get(owlsim_url, params)
        sim_resp = sim_req.json()

        disease_rank = None
        disease = disease_synth_map[synth_patient]
        if 'matches' not in sim_resp:
            logger.info(sim_profiles[synth_patient])
            continue

        result_count = len(sim_resp['matches'])
        # result_count = 7398

        # find where the disease is ranked
        for match in sim_resp['matches']:
            if match['matchId'] == disease:
                disease_rank = match['rank']

        # list to store adjusted ranks
        ranks = []

        # Resolve ties, take the average rank
        last_rank = 0
        tie_count = 0
        for match in sim_resp['matches']:
            if last_rank < match['rank']:
                if tie_count == 0:
                    ranks.append(last_rank)
                else:
                    avg_rank = round(( 2*last_rank + tie_count ) / 2)
                    tied_ranks = [avg_rank for n in range(tie_count)]
                    ranks.extend(tied_ranks)
                    # reset tie count
                    tie_count = 0
            else:
                tie_count += 1

            last_rank = match['rank']

        # make sure last
        # DRY violation, refactor
        if tie_count > 0:
            avg_rank = round((2 * last_rank + tie_count) / 2)
            tied_ranks = [avg_rank for n in range(tie_count)]
            ranks.extend(tied_ranks)

        positives = [0 for r in range(0, ranks_to_eval + 1)]
        for rank in ranks[1:]:
            positives[rank] += 1

        for rank in range(1, ranks_to_eval + 1):
            true_pos, false_pos, false_neg, true_neg = confusion_by_rank[rank]

            positive_diseases = sum(positives[0:rank+1])

            if disease_rank <= rank:
                true_pos += 1
                true_neg += result_count - positive_diseases
                false_pos += positive_diseases - 1
            else:
                false_neg += 1
                false_pos += positive_diseases
                true_neg += result_count - positive_diseases - 1

            confusion_by_rank[rank] = [true_pos, false_pos, false_neg, true_neg]

    queue.put(confusion_by_rank)

owlsim_match = "http://localhost:9000/api/match/naive-bayes-fixed-weight-two-state"
simulated_profiles: Dict[str, Set[str]] = {}
synth_to_disease: Dict[str, str] = {}

# Confusion matrix per rank
confusion_by_rank: Dict[int, List[int]] = {}

ranks_to_eval = 7398

for rank in range(1, ranks_to_eval+1):
    # tp, fp, fn, tn
    confusion_by_rank[rank] = [0, 0, 0, 0]

with gzip.open(args.simulated, 'rb') as synth_profiles:
    for line in synth_profiles:
        line = line.decode()
        if line.startswith('#'): continue
        patient, phenotype, disease = line.rstrip("\n").split("\t")[0:3]
        synth_to_disease[patient] = disease
        try:
            simulated_profiles[patient].add(phenotype)
        except KeyError:
            simulated_profiles[patient] = {phenotype}

output = open(args.output, 'w')

procs = []
queue = Queue()
result_list = []

# Split into chunks depending on args.processes
for chunk in [list(simulated_profiles.keys())[i::args.processes]
              for i in range(args.processes)]:
    proc = Process(target=create_confusion_matrix_per_rank,
                   args=(chunk, owlsim_match, simulated_profiles,
                         synth_to_disease, ranks_to_eval, queue))
    proc.start()
    procs.append(proc)

for i in range(args.processes):
    result_list.append(queue.get())

for proc in procs:
    proc.join()

for matrix in result_list:
    for rank, table in matrix.items():
        confusion_by_rank[rank][0] += table[0]
        confusion_by_rank[rank][1] += table[1]
        confusion_by_rank[rank][2] += table[2]
        confusion_by_rank[rank][3] += table[3]

for rank, confusion_matrix in confusion_by_rank.items():
    output.write("{}\t{}\t{}\t{}\t{}\n".format(
        rank,
        confusion_matrix[0],
        confusion_matrix[1],
        confusion_matrix[2],
        confusion_matrix[3]
    ))
