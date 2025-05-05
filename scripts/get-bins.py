"""
A subsumption analysis that produces a matrix output with bins of 0-5%, 6-10%, etc.
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
parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./enrichment.tsv")
args = parser.parse_args()

# i/o
output = open(args.output, "w")

lay_phenotypes = set()
gc_phenotypes = set()
lay_and_gc_phenotypes = set()
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

with open(args.mondo_labels, 'r') as mondo_labels:
    for line in mondo_labels:
        disease_id, disease_label = line.rstrip("\n").split("\t")
        mondo_disease_labels[disease_id] = disease_label

lay_and_gc_phenotypes = gc_phenotypes.union(lay_phenotypes)

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

# filter out diseases without phenotypes
associations = {d: phenos for d, phenos in associations.items() if len(phenos) > 0}

# Define bins in increments of 5%: 0-5, 6-10, 11-15, ..., 96-100
bins = [(i, i + 5) for i in range(0, 100, 5)]
bin_labels = [f"{start}-{end}" for start, end in bins]

# Initialize data structure to store binned counts
# For each bin, store counts for each category (lay, gc, lay_and_gc)
binned_counts = {
    'lay': {f"{start}-{end}": 0 for start, end in bins},
    'gc': {f"{start}-{end}": 0 for start, end in bins},
    'lay_and_gc': {f"{start}-{end}": 0 for start, end in bins}
}

counter = 0

# For each disease, calculate proportion_subsumed and add to appropriate bin
for disease, association in associations.items():
    if counter % 100 == 0:
        print(f"Processing disease {counter} of {len(associations)}")

    # Calculate proportion for each category using hpo_sim.proportion_subset
    # and place in appropriate bin
    if lay_phenotypes:  # Make sure there are lay phenotypes to avoid division by zero
        # Calculate proportion of lay phenotypes subsumed by this disease
        lay_proportion = hpo_sim.proportion_subset(association, lay_phenotypes, False) * 100
        print(lay_proportion)
        # Find which bin this falls into
        for (start, end) in bins:
            if start <= lay_proportion < end:  # Add small epsilon to include exact values
                binned_counts['lay'][f"{start}-{end}"] += 1
                break
        # If proportion is 0, add to the 0-5 bin
        if lay_proportion == 0:
            binned_counts['lay']["0-5"] += 1

    if gc_phenotypes:
        gc_proportion = hpo_sim.proportion_subset(association, gc_phenotypes, False) * 100
        for (start, end) in bins:
            if start <= gc_proportion < end:
                binned_counts['gc'][f"{start}-{end}"] += 1
                break
        if gc_proportion == 0:
            binned_counts['gc']["0-5"] += 1

    if lay_and_gc_phenotypes:
        lay_and_gc_proportion = hpo_sim.proportion_subset(association, lay_and_gc_phenotypes, False) * 100
        for (start, end) in bins:
            if start <= lay_and_gc_proportion < end:
                binned_counts['lay_and_gc'][f"{start}-{end}"] += 1
                break
        if lay_and_gc_proportion == 0:
            binned_counts['lay_and_gc']["0-5"] += 1

    counter += 1

print(f"Processed {counter} diseases")

# Write output with bins as columns and categories as rows
output.write("Category\t" + "\t".join(bin_labels) + "\n")

for category in ['lay', 'gc', 'lay_and_gc']:
    row_values = [str(binned_counts[category][bin_label]) for bin_label in bin_labels]
    output.write(f"{category}\t" + "\t".join(row_values) + "\n")

output.close()
print("Output written to", args.output)