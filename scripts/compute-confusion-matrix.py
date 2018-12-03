"""Given a file containing synthetic patient data, generates
confusion matrix per rank for the owlsim3 bayes classifier
"""
import argparse
import logging
from typing import Dict, Set, List
import gzip
import multiprocessing
from multiprocessing import Process, Queue
from phenom.utils.simulate import process_confusion_matrix_per_threshold
from phenom.model.synthetic import SyntheticProfile
import numpy

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description='Create a table of true positive/neg and false pos/neg '
                'for each rank in a binary classifier')
parser.add_argument('--simulated', '-s', type=str, required=True)
parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./confusion-matrix.tsv")
parser.add_argument('--threshold', '-t', type=str, required=False,
                    help='type of threshold; rank or probability',
                    default="probability")
parser.add_argument('--processes', '-p', type=int, required=False,
                    default=int(multiprocessing.cpu_count()/2))

args = parser.parse_args()


owlsim_match = "http://localhost:9000/api/match/naive-bayes-fixed-weight-two-state"

# Dictionaries used for constructing synthetic patient objects
simulated_profiles: Dict[str, Set[str]] = {}
synth_to_disease: Dict[str, str] = {}
synthetic_profiles: List[SyntheticProfile] = []

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

for synth_id, profile in simulated_profiles.items():
    syn_profile = SyntheticProfile(
        id=synth_id,
        phenotypes=profile,
        disease=synth_to_disease[synth_id]
    )
    synthetic_profiles.append(syn_profile)

classes_to_eval = 7344

threshold_typ = args.threshold

if threshold_typ == 'rank':
    cutoffs = range(1, classes_to_eval+1)
elif threshold_typ == 'probability':
    cutoffs = numpy.geomspace(1e-300, 1, classes_to_eval*5)
    cutoffs = numpy.flip(cutoffs)
    cutoffs = numpy.append(cutoffs, 0)
else:
    raise ValueError("{} not valid threshold".format(threshold_typ))

for cutoff in cutoffs:
    confusion_by_rank[cutoff] = [0, 0, 0, 0]

output = open(args.output, 'w')

procs = []
queue = Queue()
result_list = []

# Split into chunks depending on args.processes
for chunk in [synthetic_profiles[i::args.processes] for i in range(args.processes)]:
    proc = Process(target=process_confusion_matrix_per_threshold,
                   args=(chunk, owlsim_match, classes_to_eval,
                         cutoffs, threshold_typ, queue))
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
