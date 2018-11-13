import requests
import logging
import json
from json import JSONDecodeError
import copy
from typing import Dict, Tuple, Set, Optional, List

logger = logging.getLogger(__name__)

"""
A set of functions that interact with various Monarch APIs
such as scigraph, owlsim, solr, and the monarch app
"""

# Globals and Constants
SCIGRAPH_URL  = 'https://scigraph-ontology-dev.monarchinitiative.org/scigraph'
OWLSIM_URL    = 'https://monarchinitiative.org/owlsim/'
MONARCH_SCORE = 'https://monarchinitiative.org/score'
MONARCH_ASSOC = 'https://solr.monarchinitiative.org/solr/golr/select'
MONARCH_SEARCH = 'https://solr.monarchinitiative.org/solr/search/select'
# https://github.com/monarch-initiative/hpo-plain-index
HPO_SOLR = 'https://solr.monarchinitiative.org/solr/hpo-pl/select'

session = requests.Session()
adapter = requests.adapters.HTTPAdapter(max_retries=10)
session.mount('https://', adapter)


def get_solr_results(solr, params):
    solr_params = copy.deepcopy(params)  # don't mutate input
    result_count = solr_params['rows']
    if 'start' not in solr_params: solr_params['start'] = 0
    while solr_params['start'] < result_count:
        solr_request = session.get(solr, params=solr_params)
        response = solr_request.json()
        result_count = response['response']['numFound']
        solr_params['start'] += solr_params['rows']
        for doc in response['response']['docs']:
            yield doc


def get_clique_leader(id):
    url = SCIGRAPH_URL + '/dynamic/cliqueLeader/{}.json'.format(id)
    sci_request = session.get(url)
    response = sci_request.json()
    try:
        leader = {
            'id': response['nodes'][0]['id'],
            'label': response['nodes'][0]['lbl']
        }
    except IndexError:
        leader = {
            'id': id,
            'label': ''
        }
    return leader


def get_label(id):
    url = SCIGRAPH_URL + '/graph/{}.json'.format(id)
    sci_request = session.get(url)
    response = sci_request.json()
    try:
       label = response['nodes'][0]['lbl']
    except IndexError:
        label = ''
    return label


def get_direct_phenotypes(entity):
    phenotype_list = set()
    disease_label = ""
    params = {
        'wt': 'json',
        'rows': 200,
        'start': 0,
        'q': '*:*',
        'fl': 'subject_label, object',
        'fq': ['subject:"{0}"'.format(entity),
               'object_category:"phenotype"']
    }

    for doc in get_solr_results(MONARCH_ASSOC, params):
        phenotype_list.add(doc['object'])
        disease_label = doc['subject_label']

    return phenotype_list, disease_label


def get_annotation_sufficiency_score(id_list):
    phenotype_dictionary = dict()
    score = dict()

    phenotype_dictionary["features"] = list()
    for hp_id in id_list:
        phenotype_dictionary["features"].append({
            "id": hp_id,
            "label": "",
            "observed": "positive",
            "isPresent": "true"
        })
    phenotypes = json.dumps(phenotype_dictionary)

    params = {
        'annotation_profile': phenotypes
    }

    score_request = requests.post(MONARCH_SCORE, data=params)

    response = score_request.json()
    score['simple_score'] = response['simple_score']
    score['scaled_score'] = response['scaled_score']
    score['categorical_score'] = response['categorical_score']

    return score


def owlsim_compare(profile_a, profile_b):
    compare_url = OWLSIM_URL + "compareAttributeSets"
    params = {
        'a': profile_a,
        'b': profile_b
    }
    sim_req = session.post(compare_url, data=params)
    try:
        owlsim_results = sim_req.json()
        sim_score = owlsim_results['results'][0]['combinedScore']
    except JSONDecodeError:
        sim_score = 0
    except IndexError:
        sim_score = 0
        
    return sim_score


def get_mondo_classes() -> Dict[str, str]:
    mondo_diseases = dict()
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
    return mondo_diseases


def get_layslim() -> Tuple[Set, Set]:
    """
    Get the set of hpo classes with a lay synonym and the classes without
    """
    all_phenotypes = set()
    lay_phenotypes = set()  # terms w lay syn
    not_lay_phenotypes = set()  # terms w/o lay syn
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
    return lay_phenotypes, not_lay_phenotypes


def get_diseases_with_pheno_annotations(
        disease_classes:Optional[Dict] = None) -> Dict[str, str]:
    """
    Get a dictionary of diseases with one or more phenotype annotations
    """
    mondo_diseases = dict()  # Dict[str,str] - id: label
    facet_limit = 100000
    d2p_params = {
        'rows': 0,
        'q': '*:*',
        'fq': [
            'subject_category:disease',
            'object_category:phenotype',
            'relation:"RO:0002200"'
        ],
        'wt': 'json',
        'facet': 'true',
        'json.nl': 'arrarr',
        'facet.mincount': 1,
        'facet.limit': facet_limit,
        'facet.field': 'subject_closure'
    }
    solr_req = requests.get(MONARCH_ASSOC, params=d2p_params)
    facets = solr_req.json()
    facet_list = facets['facet_counts']['facet_fields']['subject_closure']

    if len(facet_list) > facet_limit:
        raise ValueError("Did not collect all diseases, increase facet limit")

    for field in facet_list:
        if disease_classes is not None \
                and field[0].startswith('MONDO') and field[0] in disease_classes:
            mondo_diseases[field[0]] = disease_classes[field[0]]
        elif field[0].startswith('MONDO'):
            mondo_diseases[field[0]] = ""

    return mondo_diseases


def get_disease_to_phenotype(
        include_inferred: Optional[bool] = True) -> List[Tuple[Set, Set]]:
    """
    Get disease to phenotype annotations
    Returns a list of of tuples [({d1,d2,d3}, {p1,p2,p3}),...]
    where d1,d2,d3 have phenotypes p1,p2,p3

    if include_inferred is set to False, sets will contain
    one disease and phenotype each

    :param include_inferred: if True include inferred associations
    """
    associations = []
    if include_inferred is True:
        fields = 'subject_closure,object_closure'
    else:
        fields = 'subject_closure,object'
    d2p_params = {
        'rows': 1000,
        'fl': fields,
        'q': '*:*',
        'fq': [
            'subject_category:disease',
            'object_category:phenotype',
            'relation:"RO:0002200"'
        ],
        'wt': 'json'
    }
    for assoc in get_solr_results(MONARCH_ASSOC, d2p_params):
        diseases = {dis for dis in assoc['subject_closure']
                    if dis.startswith('MONDO')}
        if include_inferred:
            phenotypes = {phe for phe in assoc['object_closure']
                          if phe.startswith('HP')}
        else:
            if not assoc['object'].startswith('HP'):
                continue
            phenotypes = {assoc['object']}
        associations.append((diseases, phenotypes))

    return associations
