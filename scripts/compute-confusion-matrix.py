"""Given a file containing synthetic patient data, generates
confusion matrix per rank for the owlsim3 bayes classifier
"""
import argparse
import logging
from typing import Dict, Set, List
import gzip
import multiprocessing
from multiprocessing import Process, Queue
from phenom.utils.simulate import rerank_by_average
import requests
import numpy

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


def create_confusion_matrix_per_threshold(
        synthetic_patients: List,
        owlsim_url: str,
        sim_profiles: Dict[str, Set[str]],
        disease_synth_map: Dict[str, str],
        num_classes: int,
        threshold_type: str,
        queue: Queue):

    confusion_by_rank: Dict[int, List[int]] = {}

    counter = 0
    total = len(synthetic_patients)

    if threshold_type == 'rank':
        for rank in range(1, classes_to_eval + 1):
            # tp, fp, fn, tn
            confusion_by_rank[rank] = [0, 0, 0, 0]
    elif threshold_type == 'probability':
        cutoffs = numpy.geomspace(1e-300, 1, classes_to_eval*5)
        cutoffs = numpy.flip(cutoffs)
        cutoffs = numpy.append(cutoffs, 0)
        for cutoff in cutoffs:
            confusion_by_rank[cutoff] = [0, 0, 0, 0]
    else:
        raise ValueError("{} not valid threshold".format(threshold_type))

    for synth_patient in synthetic_patients:

        if counter % 1000 == 0:
            logging.info("processed {} patients out of {}".format(counter, total))
        counter += 1

        # Useful for testing
        #if counter == 100:
        #    break

        params = {
            'id': sim_profiles[synth_patient],
            'limit': 8000
        }
        sim_req = requests.get(owlsim_url, params)
        sim_resp = sim_req.json()

        disease = disease_synth_map[synth_patient]
        if 'matches' not in sim_resp:
            raise ValueError("Could not find match for {}".format(sim_profiles[synth_patient]))

        # could dynamically generated num_classes here
        # num_classes = len(sim_resp['matches'])

        # store the disease rank or score
        disease_score = 0
        for disease_index, match in enumerate(sim_resp['matches']):
            if match['matchId'] == disease:
                if threshold_type == 'rank':
                    disease_score = disease_index
                elif threshold_type == 'probability':
                    disease_score = match['rawScore']

        positives = []
        cutoffs = []

        if threshold_type == 'rank':
            positives = [0 for r in range(0, num_classes + 1)]
            owlsim_ranks = [match['rank'] for match in sim_resp['matches']]
            avg_ranks = rerank_by_average(owlsim_ranks)
            disease_score = avg_ranks[disease_score]

            for rnk in avg_ranks:
                positives[rnk] += 1
            positives.pop(0)  # no ranks of 0, 1st rank is 0th index
            cutoffs = range(1, classes_to_eval + 1)

        elif threshold_type == 'probability':
            scores = numpy.array([match['rawScore'] for match in sim_resp['matches']])
            cutoffs = numpy.geomspace(1e-300, 1, classes_to_eval*5)
            cutoffs = numpy.flip(cutoffs)
            cutoffs = numpy.append(cutoffs, 0)
            positives = [0 for r in range(0, len(cutoffs))]
            for index, prob in enumerate(cutoffs):
                score_count = 0
                if index == 0:
                    positives[index] = 0
                else:
                    positives[index] = positives[index-1]
                for score in scores:
                    if score >= prob:
                        score_count += 1
                    else:
                        break
                positives[index] += score_count
                scores = scores[score_count:]

        else:
            raise ValueError("{} not valid threshold".format(threshold_type))

        for index, cutoff in enumerate(cutoffs):
            true_pos, false_pos, false_neg, true_neg = confusion_by_rank[cutoff]

            is_true_pos = False
            if threshold_type == 'rank':
                is_true_pos = disease_score <= cutoff
                positive_diseases = sum(positives[0:cutoff])
            elif threshold_type == 'probability':
                is_true_pos = disease_score >= cutoff
                positive_diseases = positives[index]
            else:
                raise ValueError

            if is_true_pos:
                true_pos += 1
                true_neg += num_classes - positive_diseases
                false_pos += positive_diseases - 1
            else:
                false_neg += 1
                false_pos += positive_diseases
                true_neg += num_classes - positive_diseases - 1

            confusion_by_rank[cutoff] = [true_pos, false_pos, false_neg, true_neg]

    queue.put(confusion_by_rank)


owlsim_match = "http://localhost:9000/api/match/naive-bayes-fixed-weight-two-state"
simulated_profiles: Dict[str, Set[str]] = {}
synth_to_disease: Dict[str, str] = {}

# Confusion matrix per rank
confusion_by_rank: Dict[int, List[int]] = {}

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

# len(list(
classes_to_eval = 7344

threshold = 'probability'

if threshold == 'rank':
    for rank in range(1, classes_to_eval+1):
        # tp, fp, fn, tn
        confusion_by_rank[rank] = [0, 0, 0, 0]
elif threshold == 'probability':
    cutoffs = numpy.geomspace(1e-300, 1, classes_to_eval*5)
    cutoffs = numpy.flip(cutoffs)
    cutoffs = numpy.append(cutoffs, 0)
    for cutoff in cutoffs:
        confusion_by_rank[cutoff] = [0, 0, 0, 0]
else:
    raise ValueError

output = open(args.output, 'w')

procs = []
queue = Queue()
result_list = []

# Split into chunks depending on args.processes
for chunk in [list(simulated_profiles.keys())[i::args.processes]
              for i in range(args.processes)]:
    proc = Process(target=create_confusion_matrix_per_threshold,
                   args=(chunk, owlsim_match, simulated_profiles,
                         synth_to_disease, classes_to_eval, threshold, queue))
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
