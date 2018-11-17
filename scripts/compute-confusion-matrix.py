"""Given a file containing synthetic patient data, generates
true positive rate and false positive rate for the owlsim3
bayes classifier
"""
import argparse
import logging
from typing import Dict, Set, List
import gzip
import requests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description='Create a table of true positive/neg and false pos/neg '
                'for each rank in a binary classifier')
parser.add_argument('--simulated', '-s', type=str, required=True)
parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./confusion-matrix.tsv")

args = parser.parse_args()

owlsim_match = "http://localhost:9000/api/match/naive-bayes-fixed-weight-two-state"
simulated_profiles: Dict[str, Set[str]] = {}
synth_to_disease: Dict[str, str] = {}

# Confusion matrix per rank
confusion_by_rank: Dict[int, List[int]] = {}

ranks_to_eval = 20

for rank in range(1, ranks_to_eval+1):
    # tp, fp, fn, tn
    confusion_by_rank[rank] = [0,0,0,0]

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
total = len(simulated_profiles.items())
counter = 0

for synth_patient, profile in simulated_profiles.items():

    if counter % 1000 == 0:
        logging.info("processed {} patients out of {}".format(counter, total))
    counter += 1

    params = {
        'id': profile,
        'limit': 300
    }
    sim_req = requests.get(owlsim_match, params)
    sim_resp = sim_req.json()

    disease_rank = None
    disease = synth_to_disease[synth_patient]
    result_count = len(sim_resp['matches'])

    # find where the disease is ranked
    for match in sim_resp['matches']:
        if match['matchId'] == disease:
            disease_rank = match['rank']  # TODO ties?
        else:
            # We could alternatively set the limit very high to get the real rank
            disease_rank = 9000

    #if disease_rank is None:
    #    logger.warn("Cannot find disease in results")
    #    disease_rank = 8000

    for rank in range(1, ranks_to_eval+1):
        true_pos, false_pos, false_neg, true_neg = confusion_by_rank[rank]
        if disease_rank <= rank:
            true_pos += 1
            true_neg += result_count - rank
            false_pos += rank - 1
        else:
            false_neg += 1
            false_pos += rank
            true_neg += result_count - rank - 1


for rank, confusion_matrix in confusion_by_rank.items():
    output.write("{}\t{}\t{}\t{}\t{}\n".format(
        rank,
        confusion_matrix[0],
        confusion_matrix[1],
        confusion_matrix[2],
        confusion_matrix[3]
    ))
