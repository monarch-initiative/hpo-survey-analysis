"""Perform enrichment of disease groups based on phenotype annotations

https://en.wikipedia.org/wiki/Fisher's_exact_test

Construct 2x2 contingency table and run fisher exact

               Mondo:123    Not Mondo:123    RowTotal
               ---------    -------------    ---
 # lay phenos  [a,          b]               sample_size
 # not lay     [c,          d]               bg_size - sample_size
               ---          ---
               nCls         nNotCls

https://github.com/biolink/ontobio/blob/d2ab6f/ontobio/assocmodel.py#L432
#from phenom.math.enrichment import fisher_exact
scipy version of fisher exact is much faster
"""
from phenom import monarch
from phenom.utils.owl_utils import get_closure, get_descendants
from scipy.stats import fisher_exact
from rdflib import Graph
import gzip
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
parser.add_argument('--output', '-o', type=str, required=False,
                    help='Location of output file', default="./enrichment.tsv")
args = parser.parse_args()

# i/o
output = open(args.output, "w")

# for including inferred phenotype annotations
# inferred disease annotations always included
include_inferred = True
lay_phenotypes = set()
not_lay_phenotypes = set()
mondo_diseases_tmp = dict()
mondo_diseases = dict()
associations = []

if args.lay_pheno:
    with open(args.lay_pheno, 'r') as lay_pheno:
        for line in lay_pheno:
            lay_phenotypes.add(line.rstrip("\n"))

else:
    lay_phenotypes, not_lay_phenotypes = monarch.get_layslim()

# all_phenotypes = lay_phenotypes | not_lay_phenotypes

if args.mondo_labels:
    with open(args.mondo_labels, 'r') as mondo_labels:
        for line in mondo_labels:
            disease_id, disease_label = line.rstrip("\n").split("\t")
            mondo_diseases_tmp[disease_id] = disease_label
else:
    mondo_diseases_tmp = monarch.get_mondo_classes()


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

    mondo = Graph()
    mondo.parse(gzip.open("../data/owl/mondo.owl.gz", 'rb'), format='xml')

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

            disease_closure = get_closure(mondo, disease, root='MONDO:0000001')

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

    if counter % 1000 == 0:
        print("Processed {} out of {} diseases".format(counter, hypothesis_count))

    # phenotypes annotated to disease class
    lay_annotated = set()
    background_annot = set()
    for assoc in associations:
        if disease in assoc[0]:
            lay_in_disease = lay_phenotypes & assoc[1]
            not_lay_dis = not_lay_phenotypes & assoc[1]
            lay_annotated = lay_annotated | lay_in_disease
            background_annot = background_annot | not_lay_dis

    # number of lay phenotypes not annotated to dis class
    lay_not_annotated = len(lay_phenotypes) - len(lay_annotated)

    # number of non-lay annotated not annotated to dis class
    background_not_annot = len(not_lay_phenotypes) - len(background_annot)

    # Run fisher exact
    # swap rows to test enrichment on terms w/o lay syn
    matrix = [
        [len(lay_annotated), lay_not_annotated],
        [len(background_annot), background_not_annot]
    ]

    odds_ratio, p_value = fisher_exact(matrix, 'greater')
    results.append({
        'id': disease,
        'label': label,
        'p-value': p_value,
        'p-value-correct': p_value * hypothesis_count
    })
    counter += 1

print("Processed {} out of {} diseases".format(counter, hypothesis_count))

associations = None
# sort results
sorted_results = sorted(results, key=lambda k: k['p-value'])

for res in sorted_results:
    output.write("{}\t{}\t{}\t{}\n".format(
        res['id'],
        res['label'],
        res['p-value'],
        res['p-value-correct']
    ))
