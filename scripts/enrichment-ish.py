
from phenom import monarch
from phenom.utils.owl_utils import get_closure, get_descendants
from rdflib import Graph
import argparse

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

# for including inferred phenotype annotations
# inferred disease annotations always included
include_inferred = True
lay_phenotypes = set()
gc_phenotypes = set()
mondo_leaf_nodes = set()
not_lay_phenotypes = set()
mondo_diseases_tmp = dict()
mondo_diseases = dict()
associations = []

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
        mondo_diseases_tmp[disease_id] = disease_label

hpo = Graph()
hpo.parse("../data/owl/hp.owl", format='xml')
root = "HP:0000118"
phenotype_terms = get_descendants(hpo, root)
not_lay_phenotypes = phenotype_terms - lay_phenotypes

if not include_inferred:
    hpo = None  # free up some mem

print("{} phenotypes in outer set".format(len(not_lay_phenotypes)))

if args.mondo_assoc:
    print("Fetching associations from cache file")
    counter = 0

    #mondo = Graph()
    #mondo.parse(gzip.open("../data/owl/mondo.owl.gz", 'rb'), format='xml')

    with open(args.mondo_assoc, 'r') as mondo_labels:
        for line in mondo_labels:
            if line.startswith('#'): continue
            if not line.startswith('MONDO'): continue
            if counter % 10000 == 0:
                print("Processed {} associations".format(counter))
            disease, phenotype = line.rstrip("\n").split("\t")[0:2]
            try:
                mondo_diseases[disease] = mondo_diseases_tmp[disease]
            except KeyError:
                mondo_diseases[disease] = "obsoleted class"

            disease_closure = {disease}

            for dis in disease_closure:
                try:
                    mondo_diseases[dis] = mondo_diseases_tmp[dis]
                except KeyError:
                    mondo_diseases[dis] = "obsoleted class"

            if include_inferred:
                phenotype_closure = get_closure(hpo, phenotype, root='HP:0000118')
                associations.append((disease_closure, phenotype_closure))
            else:
                associations.append((disease_closure, {phenotype}))
            counter += 1

else:
    print("Fetching associations from golr")
    # Get associations using golr
    mondo_diseases = monarch.get_diseases_with_pheno_annotations(mondo_diseases_tmp)
    associations = monarch.get_disease_to_phenotype(include_inferred)

print("Finished fetching associations")
mondo_diseases_tmp = None

# number of hypotheses for bonferroni correction
# hypotheses = number of disease classes with at least 1 association
hypothesis_count = len(mondo_diseases.keys())

results = []
counter = 0

for disease, label in mondo_diseases.items():

    if not disease in mondo_leaf_nodes: continue

    if counter % 10 == 0:
        print("Processed {} out of {} diseases".format(counter, hypothesis_count))

    # phenotypes annotated to disease class
    lay_annotated = set()
    lay_other = set()
    gc_other = set()
    gc_annotated = set()
    all_annot = set()
    for assoc in associations:
        if disease in assoc[0]:
            lay_in_disease = lay_phenotypes & assoc[1]
            gc_in_disease = gc_phenotypes & assoc[1]
            lay_annotated = lay_annotated | lay_in_disease
            gc_annotated = gc_annotated | gc_in_disease
            all_annot = all_annot | assoc[1]
        else:
            lay_in_other = lay_phenotypes & assoc[1]
            gc_in_other = gc_phenotypes & assoc[1]
            lay_other = lay_in_other | lay_other
            gc_other = gc_in_other | gc_other

    results.append({
        'id': disease,
        'label': label,
        'clinical annotated to disease': len(all_annot),
        'lay annotated to disease': len(lay_annotated),
        'gc annotated to disease': len(gc_annotated),
        'lay other': len(lay_other),
        'gc other': len(gc_other)
    })
    counter += 1

print("Processed {} out of {} diseases".format(counter, hypothesis_count))

associations = None
# sort results
# sorted_results = sorted(results, key=lambda k: k['p-value'])

headers = [
    'id',
    'label',
    'clinical annotated to disease',
    'lay annotated to disease',
    'gc annotated to disease',
    'lay other',
    'gc other'
]

header_line = "\t".join(headers)

output.write(f"{header_line}\n")

for res in results:
    result = "\t".join([
        str(res[field]) for field in headers
    ])
    output.write(f"{result}\n")
