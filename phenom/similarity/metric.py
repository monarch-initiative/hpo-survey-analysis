from typing import Set, Union, Iterable, Dict, Optional
from itertools import chain
from phenom.utils import owl_utils
from phenom.math import math
from rdflib import RDFS, Graph, URIRef


# Union types
Num = Union[int, float]


def jaccard(set1: Set, set2: Set) -> float:
    return len(set1.intersection(set2))/len(set1.union(set2))


def pairwise_jaccard(pheno_a: str, pheno_b: str, graph: Graph, root:str) -> float:
    predicate = RDFS['subClassOf']
    return jaccard(
        owl_utils.get_closure(graph, pheno_a, predicate, root),
        owl_utils.get_closure(graph, pheno_b, predicate, root)
    )


def get_mica_ic(
        pheno_a: str,
        pheno_b: str,
        graph: Graph,
        ic_map: Dict[str, float],
        root) -> float:
    predicate = RDFS['subClassOf']
    p1_closure = owl_utils.get_closure(graph, pheno_a, predicate, root)
    p2_closure = owl_utils.get_closure(graph, pheno_b, predicate, root)
    return max([ic_map[parent]for parent in p1_closure.intersection(p2_closure)])


def jac_ic_geomean(
        pheno_a: str,
        pheno_b: str,
        graph: Graph,
        ic_map: Dict[str, float],
        root) -> float:
    jaccard_sim = pairwise_jaccard(pheno_a, pheno_b, graph, root)
    mica = get_mica_ic(pheno_a, pheno_b, graph, ic_map, root)
    return math.geometric_mean([jaccard_sim, mica])


def get_profile_closure(
        profile: Iterable[str],
        graph: Graph,
        root: str,
        predicate: Optional[URIRef] = RDFS['subClassOf'],
        negative: Optional[bool] = False) -> Set[str]:
    """
    Given a list of phenotypes, get the reflexive closure for each phenotype
    stored in a single set.  This can be used for jaccard similarity or
    simGIC
    """
    return set(chain.from_iterable(
        [owl_utils.get_closure(
            graph, pheno, predicate, root, reflexive=True, negative=negative)
         for pheno in profile])
    )
