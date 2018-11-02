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
"""
from phenom.monarch import get_solr_results
from phenom.math.enrichment import fisher_exact
#from scipy.stats import fisher_exact

# monarch associations
MONARCH_ASSOC = 'https://solr.monarchinitiative.org/solr/golr/select'

# monarch search
MONARCH_SEARCH = 'https://solr.monarchinitiative.org/solr/search/select'

# https://github.com/monarch-initiative/hpo-plain-index
HPO_SOLR = 'https://solr.monarchinitiative.org/solr/hpo-pl/select'

# i/o
output = open("/path/to/lay-enrich.tsv", "w")

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

# Get all diseases using monarch search index
# ?q=*:*&fl=id,label&wt=json&fq=prefix:MONDO
params = {
    'rows': 100,
    'fl': 'id,label',
    'q': '*:*',
    'fq': 'prefix:MONDO',
    'wt': 'json'
}
for disease in get_solr_results(MONARCH_SEARCH, params):
    mondo_diseases[disease['id']] = disease['label'][0]

# Get associations using golr
params = {
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
associations = []
for assoc in get_solr_results(MONARCH_ASSOC, params):
    disease_set = {dis for dis in assoc['subject_closure']
                   if dis.startswith('MONDO')}
    phenotype_set = {phe for phe in assoc['object_closure']
                     if phe.startswith('HP')}
    #if not assoc['object'].startswith('HP'): continue
    #phenotype_set = {assoc['object']}
    associations.append((disease_set, phenotype_set))

# number of hypotheses for bonferroni correction
# hypotheses = number of disease classes
hypothesis_count = len(mondo_diseases.keys())

# background count = number of phenotype classes
background_count = len(all_phenotypes)

results = []

for disease, label in mondo_diseases.items():

    # phenotypes annotated to disease class
    lay_annotated = set()
    background_annot = set()
    for assoc in associations:
        if disease in assoc[0]:
            lay_in_disease = lay_phenotypes & assoc[1]
            not_lay_dis = not_lay_phenotypes & assoc[1]
            if len(lay_in_disease) > 0:
                lay_annotated = lay_annotated | lay_in_disease
            elif len(not_lay_dis) > 0:
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

    p_value = fisher_exact(matrix, 'greater')
    results.append({
        'id': disease,
        'label': label,
        'p-value': p_value,
        'p-value-correct': p_value * hypothesis_count
    })

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
