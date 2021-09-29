"""
A subsumption analysis that includes more data points than subsumption-counts.py
"""
import argparse
from collections import defaultdict
from pumpkin_py import build_graph_from_rdflib, GraphSemSim

parser = argparse.ArgumentParser(
    description='Perform disease enrichment using fisher exact test'
)
parser.add_argument('--mondo_labels', '-m', type=str, required=False,
                    help='Cached 2 column diseases and labels tsv')
parser.add_argument('--mondo_assoc', '-a', type=str, required=False,
                    help='path to mondo to phenotype assoc 2 column tsv')
parser.add_argument('--lay_pheno', '-l', type=str, required=False,
                    help='path to lay phenotypes 1 column txt')
parser.add_argument('--gc_pheno', '-g', type=str, required=False,
                    help='path to gc phenotypes 1 column txt')
parser.add_argument('--mondo_leaf', '-ml', type=str, required=False,
                    help='path to gc phenotypes 1 column txt')
parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./enrichment.tsv")
args = parser.parse_args()

# i/o
output = open(args.output, "w")

lay_phenotypes = set()
gc_phenotypes = set()
mondo_leaf_nodes = set()
hpo_with_annotation = defaultdict(int)
mondo_disease_labels = dict()
mondo_diseases = dict()
associations = defaultdict(set)

with open(args.lay_pheno, 'r') as lay_pheno:
    for line in lay_pheno:
        lay_phenotypes.add(line.rstrip("\n"))

with open(args.gc_pheno, 'r') as gc_pheno:
    for line in gc_pheno:
        gc_phenotypes.add(line.rstrip("\n"))

with open(args.mondo_leaf, 'r') as mondo_leaf:
    for line in mondo_leaf:
        mondo_leaf_nodes.add(line.rstrip("\n"))

with open(args.mondo_labels, 'r') as mondo_labels:
    for line in mondo_labels:
        disease_id, disease_label = line.rstrip("\n").split("\t")
        mondo_disease_labels[disease_id] = disease_label


root = "HP:0000118"
hpo = build_graph_from_rdflib("../data/owl/hp.owl", root)
hpo_sim = GraphSemSim(hpo)

counter = 0

print("Fetching associations from cache file")

with open(args.mondo_assoc, 'r') as mondo_labels:
    for line in mondo_labels:
        if line.startswith('#'): continue
        if not line.startswith('MONDO'): continue
        if counter % 10000 == 0:
            print("Processed {} associations".format(counter))
        disease, phenotype = line.rstrip("\n").split("\t")[0:2]
        try:
            mondo_diseases[disease] = mondo_disease_labels[disease]
        except KeyError:
            mondo_diseases[disease] = "obsoleted class"

        phenotype_closure = {hpo.id_map.inverse[p] for p in hpo.get_closure(phenotype)}
        if 'HP:0000118' not in phenotype_closure:
            # Filter out non phenotype terms
            print("closure without 118")
            continue
        for phen in phenotype_closure:
            hpo_with_annotation[phen] += 1

        associations[disease] = associations[disease].union(phenotype_closure)

        counter += 1


print(f"Processed {counter} associations")

results = []
counter = 0
hpo = None

hpo_with_three_or_more_annotations = {k for k,v in hpo_with_annotation.items() if v >= 3}
hpo_with_five_or_more_annotations = {k for k,v in hpo_with_annotation.items() if v >= 5}

# filter out diseases without phenotypes
associations = {disease: associations for disease, associations in associations.items() if len(associations) > 0}

for disease, association in associations.items():

    if counter % 100 == 0:
        print("Processed {} diseases".format(counter))

    # phenotypes annotated to disease class
    lay_annotated = set()
    lay_other = set()
    gc_other = set()
    gc_annotated = set()
    all_annot = set()

    lay_subsumed_gt_25 = 0
    lay_subsumed_gt_50 = 0
    lay_subsumed_gt_75 = 0
    lay_subsumed_100 = 0

    gc_subsumed_gt_25 = 0
    gc_subsumed_gt_50 = 0
    gc_subsumed_gt_75 = 0
    gc_subsumed_100 = 0

    lay_in_disease = lay_phenotypes & association
    gc_in_disease = gc_phenotypes & association
    lay_annotated = lay_annotated | lay_in_disease
    gc_annotated = gc_annotated | gc_in_disease
    all_annot = all_annot | association

    for other_disease, other_association in associations.items():
        if other_disease == disease: continue

        # len(intersection(lay, gold)) / len(lay)
        lay_proportion_subsumed = hpo_sim.proportion_subset(lay_in_disease, other_association, False)

        # len(intersection(lay, gold)) / len(gold)
        # lay_proportion_subsumed = hpo_sim.proportion_subset(other_association, lay_in_disease, False)

        # why is there no case switch
        if lay_proportion_subsumed > .25:
            lay_subsumed_gt_25 += 1

        if lay_proportion_subsumed > .50:
            lay_subsumed_gt_50 += 1

        if lay_proportion_subsumed > .75:
            lay_subsumed_gt_75 += 1

        if lay_proportion_subsumed == 1.0:
            lay_subsumed_100 += 1

        # len(intersection(gc, gold)) / len(gc)
        gc_proportion_subsumed = hpo_sim.proportion_subset(gc_in_disease, other_association, False)

        # len(intersection(gc, gold)) / len(gold)
        # gc_proportion_subsumed = hpo_sim.proportion_subset(other_association, gc_in_disease, False)

        if gc_proportion_subsumed > .25:
            gc_subsumed_gt_25 += 1

        if gc_proportion_subsumed > .50:
            gc_subsumed_gt_50 += 1

        if gc_proportion_subsumed > .75:
            gc_subsumed_gt_75 += 1

        if gc_proportion_subsumed == 1.0:
            gc_subsumed_100 += 1


    results.append({
        'id': disease,
        'label': mondo_disease_labels[disease],
        'clinical annotated to disease': len(all_annot),
        'lay annotated to disease': len(lay_annotated),
        'gc annotated to disease': len(gc_annotated),
        'lay two': len(lay_annotated & hpo_with_three_or_more_annotations),
        'gc two': len(gc_annotated & hpo_with_three_or_more_annotations),
        'lay four': len(lay_annotated & hpo_with_five_or_more_annotations),
        'gc four': len(gc_annotated & hpo_with_five_or_more_annotations),
        'lay 25': lay_subsumed_gt_25,
        'gc 25': gc_subsumed_gt_25,
        'lay 50': lay_subsumed_gt_50,
        'gc 50': gc_subsumed_gt_50,
        'lay 75': lay_subsumed_gt_75,
        'gc 75': gc_subsumed_gt_75,
        'lay 100': lay_subsumed_100,
        'gc 100': gc_subsumed_100,

    })
    counter += 1

print("Processed {} diseases".format(counter))

headers = [
    'id',
    'label',
    'clinical annotated to disease',
    'lay annotated to disease',
    'gc annotated to disease',
    'lay two',
    'gc two',
    'lay four',
    'gc four',
    'lay 25',
    'gc 25',
    'lay 50',
    'gc 50',
    'lay 75',
    'gc 75',
    'lay 100',
    'gc 100',
]

header_line = "\t".join(headers)

output.write(f"{header_line}\n")

for res in results:
    result = "\t".join([
        str(res[field]) for field in headers
    ])
    output.write(f"{result}\n")
