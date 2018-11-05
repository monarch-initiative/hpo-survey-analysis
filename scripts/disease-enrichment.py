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
from phenom.monarch import get_solr_results, get_mondo_classes
from scipy.stats import fisher_exact
import requests
from copy import deepcopy

# monarch associations
MONARCH_ASSOC = 'https://solr.monarchinitiative.org/solr/golr/select'

# https://github.com/monarch-initiative/hpo-plain-index
HPO_SOLR = 'https://solr.monarchinitiative.org/solr/hpo-pl/select'

# i/o
output = open("/path/to/out.tsv", "w")

all_phenotypes = set()
lay_phenotypes = set()  # terms w lay syn
not_lay_phenotypes = set()  # terms w/o lay syn
mondo_diseases = dict()  # Dict[str,str] - id: label

# Get all phenotypes in HPO as background using phenotypr index
params = {
    'rows': 100,
    'fl': 'id,has_pl_syn',
    'q': '*:*'
}

for pheno in get_solr_results(HPO_SOLR, params):
    all_phenotypes.add(pheno['id'])
    if pheno['has_pl_syn'] is True:
        lay_phenotypes.add(pheno['id'])

not_lay_phenotypes = all_phenotypes - lay_phenotypes

mondo_diseases_tmp = get_mondo_classes()

# Get all disease classes with direct or inferred phenotype annotations
facet_limit = 100000
d2p_params = {
    'rows': 1000,
    'fl': 'subject_closure,object_closure',  # object for direct
    'q': '*:*',
    'fq': [
        'subject_category:disease',
        'object_category:phenotype',
        'relation:"RO:0002200"'
    ],
    'wt': 'json'
}
facet_params = deepcopy(d2p_params)
facet_params['rows'] = 0
facet_params['facet'] = "true"
facet_params['json.nl'] = "arrarr"
facet_params['facet.mincount'] = "1"
facet_params['facet.limit'] = facet_limit
facet_params['facet.field'] = "subject_closure"

solr_req = requests.get(MONARCH_ASSOC, params=facet_params)
facets = solr_req.json()
facet_list = facets['facet_counts']['facet_fields']['subject_closure']

if len(facet_list) > facet_limit:
    raise ValueError("Did not collect all diseases, increase facet limit")

for field in facet_list:
    if field[0].startswith('MONDO') and field[0] in mondo_diseases_tmp:
        mondo_diseases[field[0]] = mondo_diseases_tmp[field[0]]

mondo_diseases_tmp = None

# Get associations using golr
associations = []
for assoc in get_solr_results(MONARCH_ASSOC, d2p_params):
    diseases = {dis for dis in assoc['subject_closure']
                if dis.startswith('MONDO')}
    phenotypes = {phe for phe in assoc['object_closure']
                  if phe.startswith('HP')}
    #if not assoc['object'].startswith('HP'): continue
    #phenotypes = {assoc['object']}
    associations.append((diseases, phenotypes))

# number of hypotheses for bonferroni correction
# hypotheses = number of disease classes with at least 1 association
hypothesis_count = len(mondo_diseases.keys())

# background count = number of phenotype classes
background_count = len(all_phenotypes)

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
