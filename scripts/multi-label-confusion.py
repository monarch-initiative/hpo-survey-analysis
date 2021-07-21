"""Given a file containing synthetic patient data, generates
confusion matrix per disease, treating diseases as labels

https://scikit-learn.org/stable/modules/generated/sklearn.metrics.multilabel_confusion_matrix.html
https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html
"""
import argparse
import logging
from typing import Dict, Set, List
import gzip
from phenom.model.synthetic import SyntheticProfile
import numpy
from collections import defaultdict
from pathlib import Path
from pumpkin_py import ICSemSim, flat_to_annotations, build_ic_graph_from_closures


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description='Create a table of true positive/neg and false pos/neg '
                'for each rank in a binary classifier')
parser.add_argument('--simulated', '-s', type=str, required=True)
parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./confusion-matrix.tsv")


args = parser.parse_args()

output = open(args.output, 'w')
annotation_file = Path(__file__).parents[1] /  'data' / 'synthetic' / 'gold-standard.tsv.gz'
closures = Path(__file__).parents[1] / 'data' / 'hp-closures.tsv'

root = "HP:0000118"

diseases = set()

print("Loading closures")

with gzip.open(annotation_file, 'rt') as annot_file:
    annotations = flat_to_annotations(annot_file)

with gzip.open(annotation_file, 'rt') as annot_file:
    for line in annot_file:
        row = line.rstrip().split('\t')
        diseases.add(row[0])

print("Building graph")
with open(closures, 'r') as closure_file:
    graph = build_ic_graph_from_closures(closure_file, root, annotations)


ic_sim = ICSemSim(graph)


# Dictionaries used for constructing synthetic patient objects
simulated_profiles: Dict[str, Set[str]] = defaultdict(set)
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
        simulated_profiles[patient].add(phenotype)

for synth_id, profile in simulated_profiles.items():
    syn_profile = SyntheticProfile(
        id=synth_id,
        phenotypes=profile,
        disease=synth_to_disease[synth_id]
    )
    synthetic_profiles.append(syn_profile)

cutoffs = numpy.arange(0, 100, 0.2)
cutoffs = numpy.flip(cutoffs)

# An object of diseases as keys
# values are tuples of threshold, confusion matrix
disease_confusion_matrices = {
    disease: [(cutoff, [0, 0, 0, 0]) for cutoff in cutoffs]
    for disease in diseases
}

print("Generating confusion matrices")
for patient_index, syn_profile in enumerate(synthetic_profiles):
    if patient_index % 1000 == 0:
        print(f"Processed {patient_index} profiles out of {len(synthetic_profiles)}")

    for test_disease, gold_standard in annotations.items():
        phenodigm_score = ic_sim.phenodigm_compare(syn_profile.phenotypes, gold_standard)
        for index, confusion in enumerate(disease_confusion_matrices[test_disease]):
            confusion_threshold, confusion_matrix = confusion
            true_pos, false_pos, false_neg, true_neg = confusion_matrix
            if test_disease == syn_profile.disease:
                # if the score is greater than the threshold this is a true positive
                # else it's false negative
                if phenodigm_score >= confusion_threshold:
                    true_pos += 1
                else:
                    false_neg += 1
            else:
                # if the score is greater than the threshold this is a false positive
                # else it's a true negative
                if phenodigm_score >= confusion_threshold:
                    false_pos += 1
                else:
                    true_neg += 1
            confusion_matrix = [true_pos, false_pos, false_neg, true_neg]
            disease_confusion_matrices[test_disease][index] = (confusion_threshold, confusion_matrix)

for disease, confusions in disease_confusion_matrices.items():
    for confusion in confusions:
        confusion_threshold, confusion_matrix = confusion
        true_pos, false_pos, false_neg, true_neg = confusion_matrix
        output_line = "\t".join([
            disease,
            str(confusion_threshold),
            str(true_pos),
            str(false_pos),
            str(false_neg),
            str(true_neg)
        ])
        output.write(f'{output_line}\n')
