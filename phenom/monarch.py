import requests
import logging
import json
from json import JSONDecodeError
import copy

logger = logging.getLogger(__name__)

"""
A set of functions that interact with various Monarch APIs
such as scigraph, owlsim, solr, and the monarch app
"""

# Globals and Constants
SCIGRAPH_URL =    'https://scigraph-ontology-dev.monarchinitiative.org/scigraph'
OWLSIM_URL =      'https://monarchinitiative.org/owlsim/'
MONARCH_SCORE =   'https://monarchinitiative.org/score'
MONARCH_ASSOC =   'https://solr.monarchinitiative.org/solr/golr/select'

session = requests.Session()
adapter = requests.adapters.HTTPAdapter(max_retries=10)
session.mount('https://', adapter)


def get_solr_results(solr, params):
    solr_params = copy.deepcopy(params)  # don't mutate input
    result_count = solr_params['rows']
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

