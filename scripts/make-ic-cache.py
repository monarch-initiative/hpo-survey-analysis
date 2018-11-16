from typing import Dict, Set, List, Tuple
import csv
import requests
from contextlib import closing
from phenom import monarch
from phenom.utils import owl_utils
from phenom.math import math_utils
from prefixcommons import contract_uri, expand_uri
from prefixcommons.curie_util import NoPrefix
from rdflib import Graph, URIRef, OWL, RDFS
import logging
import argparse
import gzip
from argparse import ArgumentError

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# globals and constants
#HPO_ONTOLOGY = "http://purl.obolibrary.org/obo/hp.owl"
HPO_ONTOLOGY = "../data/hp.owl"

# contains manual and semi-automated annotations created by the HPO-team.
# These are annotations of OMIM-, Orphanet-, and DECIPHER-entries
#HPO_DATA = "http://compbio.charite.de/jenkins/job/hpo.annotations/" \
#           "lastStableBuild/artifact/misc/phenotype_annotation.tab"

HPO_DATA = "../data/phenotype_annotation.tab"


def main():

    parser = argparse.ArgumentParser(
        description='Generate information content for each HPO class using the '
                    'HPO phenotype annotation file ')
    parser.add_argument('--mondo_cache', '-m', type=str, required=False,
                        help='Cached 2 column disease phenotype tsv')
    parser.add_argument('--mondo_output', '-mo', type=str, required=False,
                        help='path to mondo 2 column file (if no cache)')
    parser.add_argument('--output', '-o', type=str, required=False,
                        help='Location of output file', default="./ic-cache.tsv")
    args = parser.parse_args()

    if not args.mondo_cache and not args.mondo_output:
        raise ArgumentError(
            args.mondo_cache, 'One of --mondo_cache or '
                              '--mondo_output must be provided')

    # i/o
    output_file = open(args.output, 'w')

    explicit_annotations = 0
    pheno_annotations: Dict[str, int] = {}
    hp_graph = Graph()
    hp_graph.parse(HPO_ONTOLOGY, format='xml')
    root = "HP:0000118"
    all_hpo = "HP:0000001"

    phenotype_terms = owl_utils.get_descendants(hp_graph, root)
    all_terms = owl_utils.get_descendants(hp_graph, all_hpo)
    pheno_annotations = {phenotype: 0 for phenotype in phenotype_terms}

    # filter out clinical course, inheritance, modifiers, etc
    non_phenotype_terms = all_terms.difference(phenotype_terms)

    if not args.mondo_cache:
        args.mondo_cache = args.mondo_output
        mondo_cache = open(args.mondo_output, 'w')
        header = "#" + "\t".join([
            "disease", "phenotype", "onset", "frequency", "severity"])
        mondo_cache.write(header + "\n")
        d2p_map = _process_hpo_data(HPO_DATA)
        for d2p, pheno_info in d2p_map.items():
            disease, phenotype = d2p.split("-")
            mondo_cache.write("{}\t{}\t{}\n".format(
                disease, phenotype, "\t".join(pheno_info))
            )
        mondo_cache.close()

    with open(args.mondo_cache, 'r') as cache_file:
        reader = csv.reader(cache_file, delimiter='\t', quotechar='\"')
        for row in reader:
            if row[0].startswith('#'): continue
            (mondo_id, phenotype_id) = row[0:2]
            if phenotype_id in non_phenotype_terms:
                continue
            explicit_annotations += 1
            for pheno in owl_utils.get_ancestors(hp_graph, phenotype_id, root=root):
                pheno_annotations[pheno] += 1

    for phenotype, annot_count in pheno_annotations.items():
        output_file.write("{}\t{}\n".format(
            phenotype,
            math_utils.information_content(annot_count / explicit_annotations)
        ))


def _process_hpo_data(file_path: str) -> Dict[str, List[str]]:
    logger.info("loading mondo into memory")
    mondo = Graph()
    mondo.parse(gzip.open("../data/mondo.owl.gz", 'rb'), format='xml')
    logger.info("finished loading mondo")

    mondo_merged_lines: List[str] = []
    disease_info: Dict[str, List[str]] = {}

    if file_path.startswith("http"):
        context_manager = closing(requests.get(file_path))
    else:
        context_manager = open(file_path, "r")

    # https://stackoverflow.com/a/35371451
    with context_manager as file:
        if file_path.startswith("http"):
            file = file.content.decode('utf-8').splitlines()
        reader = csv.reader(file, delimiter='\t', quotechar='\"')
        counter = 0
        for row in reader:
            try:
                (db, num, name, severity, pheno_id, publist, eco, onset, freq) = row[0:9]
            except ValueError:
                logger.warning("Too few values in row {}".format(row))
                continue

            # Align Id prefixes
            if db == 'MIM':      db = 'OMIM'
            if db == 'ORPHA':    db = 'Orphanet'
            if db == 'ORPHANET': db = 'Orphanet'

            disease_id = "{}:{}".format(db, num)
            disease_iri = URIRef(expand_uri(disease_id, strict=True))
            mondo_curie = None
            mondo_iri = None
            for subj in mondo.subjects(OWL['equivalentClass'], disease_iri):
                curie = contract_uri(str(subj), strict=True)[0]
                if curie.startswith('MONDO'):
                    mondo_curie = curie
                    mondo_iri = subj
                    break
            if mondo_curie is None:
                logger.warn("No mondo id for {}".format(disease_id))
                continue

            has_omim = False
            for obj in mondo.objects(mondo_iri, OWL['equivalentClass']):
                try:
                    curie = contract_uri(str(obj), strict=True)[0]
                except NoPrefix:
                    continue
                if curie.startswith('OMIM'):
                    has_omim = True

            # use scigraph instead of the above
            # mondo_node = monarch.get_clique_leader(disease_id)
            # mondo_curie = mondo_node['id']
            if mondo_curie is not None and 'hgnc' in mondo_curie:
                # to keep these, likely decipher IDs
                # mondo_curie = disease_id
                continue

            if disease_id.startswith('Orphanet') \
                    and has_omim is False \
                    and len(list(mondo.objects(mondo_iri, RDFS['subClassOf']))) > 0:
                # disease is a disease group, skip
                logger.info("{} is a disease group, skipping".format(disease_id))
                continue

            mondo_merged_lines.append((mondo_curie, pheno_id, onset, freq, severity))

            counter += 1
            if counter % 10000 == 0:
                logger.info("processed {} rows".format(counter))

    logger.info("processed {} rows".format(counter))

    for line in mondo_merged_lines:
        key = "{}-{}".format(line[0], line[1])
        values = [line[2], line[3], line[4]]
        if key in disease_info and disease_info[key] != values:
            logger.warning("Metadata for {} and {} mismatch: {} vs {}".format(
                line[0], line[1], values, disease_info[key])
            )
            # attempt to merge by collapsing freq, onset, severity
            # that is empty in one disease but not another
            # conflicts will defer to the disease first inserted
            merged_disease_info = disease_info[key]
            for index, val in enumerate(values):
                if val == disease_info[key][index] \
                        or val == '' and disease_info[key][index] != '':
                    continue
                elif val != '' and disease_info[key][index] == '':
                    merged_disease_info[index] = val
                else:
                    logger.warning("Cannot merge {} and {} for {}".format(
                        values, disease_info[key], line[0])
                    )
        else:
            disease_info[key] = values

    return disease_info


if __name__ == "__main__":
    main()
