"""
A subsumption analysis that produces a matrix output with threshold bins (0-5, 6-10, etc.)
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

counter = 0

hpo_with_three_or_more_annotations = {k for k, v in hpo_with_annotation.items() if v >= 3}
hpo_with_five_or_more_annotations = {k for k, v in hpo_with_annotation.items() if v >= 5}

# filter out diseases without phenotypes
associations = {disease: phenos for disease, phenos in associations.items() if len(phenos) > 0}

# Define threshold bins
# We'll use 0-5, 6-10, 11-15, ..., 96-100
bins = [(i, i + 5) for i in range(0, 101, 5)]

# Initialize dictionaries to hold binned counts for each category across all diseases
lay_bin_counts = defaultdict(int)
gc_bin_counts = defaultdict(int)
lay_and_gc_bin_counts = defaultdict(int)

for disease, association in associations.items():
    if counter % 100 == 0:
        print("Processed {} diseases".format(counter))

    lay_in_disease = lay_phenotypes & association
    gc_in_disease = gc_phenotypes & association
    lay_and_gc_in_disease = lay_and_gc_phenotypes & association

    # Skip processing if there are no phenotype annotations
    if not lay_in_disease and not gc_in_disease and not lay_and_gc_in_disease:
        counter += 1
        continue

    for other_disease, other_association in associations.items():
        if other_disease == disease: continue

        # Calculate proportion subsumed for each category
        lay_proportion = hpo_sim.proportion_subset(lay_in_disease, other_association,
                                                   False) * 100 if lay_in_disease else 0
        gc_proportion = hpo_sim.proportion_subset(gc_in_disease, other_association, False) * 100 if gc_in_disease else 0
        lay_and_gc_proportion = hpo_sim.proportion_subset(lay_and_gc_in_disease, other_association,
                                                          False) * 100 if lay_and_gc_in_disease else 0

        # Assign to appropriate bin
        for bin_start, bin_end in bins:
            bin_label = f"{bin_start}-{bin_end}"

            if bin_start <= lay_proportion < bin_end or (bin_end == 100 and lay_proportion == 100):
                lay_bin_counts[bin_label] += 1

            if bin_start <= gc_proportion < bin_end or (bin_end == 100 and gc_proportion == 100):
                gc_bin_counts[bin_label] += 1

            if bin_start <= lay_and_gc_proportion < bin_end or (bin_end == 100 and lay_and_gc_proportion == 100):
                lay_and_gc_bin_counts[bin_label] += 1

    counter += 1

print("Processed {} diseases".format(counter))

# Generate bin labels for columns
bin_labels = [f"{bin_start}-{bin_end}" for bin_start, bin_end in bins]

# Write the headers
header_line = "Category\t" + "\t".join(bin_labels)
output.write(f"{header_line}\n")

# Write the data for each category
output.write(f"lay\t" + "\t".join([str(lay_bin_counts[bin_label]) for bin_label in bin_labels]) + "\n")
output.write(f"gc\t" + "\t".join([str(gc_bin_counts[bin_label]) for bin_label in bin_labels]) + "\n")
output.write(f"lay_and_gc\t" + "\t".join([str(lay_and_gc_bin_counts[bin_label]) for bin_label in bin_labels]) + "\n")

output.close()
print("Results written to", args.output)
