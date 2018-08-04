from typing import Set, List, Optional, Dict, Sequence
from phenom.decorators import memoized
from phenom.utils import math_utils
from rdflib import URIRef, BNode, Literal, Graph, RDFS


def get_closure(
        graph: Graph,
        node: str,
        edge: Optional[URIRef]=None,
        root: Optional[str]=None,
        reflexive: Optional[bool] = True) -> Set[str]:

    hp_id = int(node.replace("HP:", ""))
    closure = get_closure_memoized(graph, hp_id, edge, root, reflexive)
    return {"HP:" + str(node).zfill(7) for node in closure}


def _get_closure(
        graph: Graph,
        node: str,
        edge: Optional[URIRef]=None,
        root: Optional[str]=None,
        reflexive: Optional[bool] = True) -> Set[str]:

    nodes = set()
    root_seen = {}
    if node is not None:
        node = hp_curie2iri(node)
    if root is not None:
        root = hp_curie2iri(root)
        root_seen = {root: 1}
    for obj in graph.transitive_objects(node, edge, root_seen):
        if isinstance(obj, Literal) or isinstance(obj, BNode):
            continue
        if not reflexive and node == obj:
            continue
        nodes.add(hp_iri2curie(str(obj)))

    # Add root to graph
    if root is not None:
        nodes.add(hp_iri2curie(str(root)))

    return nodes


@memoized
def get_closure_memoized(
        graph: Graph,
        node: int,
        edge: Optional[URIRef]=None,
        root: Optional[str]=None,
        reflexive: Optional[bool] = True) -> List[int]:

    node = "HP:" + str(node).zfill(7)
    closure = _get_closure(graph, node, edge, root, reflexive)
    return [int(id.replace("HP:","")) for id in closure]


def get_mica_ic(
        pheno_a: str,
        pheno_b: str,
        graph: Graph,
        ic_map:Dict[str, float],
        root) -> float:
    predicate = RDFS['subClassOf']
    p1_closure = get_closure(graph, pheno_a, predicate, root)
    p2_closure = get_closure(graph, pheno_b, predicate, root)
    return max([ic_map[parent]for parent in p1_closure.intersection(p2_closure)])


def get_mica_id(
        pheno_a: str,
        pheno_b: str,
        graph: Graph,
        ic_map:Dict[str, float],
        root: str) -> float:
    """
    Return ID of most informative common anscestor of two phenotypes
    Currently does not handle ambiguity (>1 equal MICAs)
    """
    predicate = RDFS['subClassOf']
    p1_closure = get_closure(graph, pheno_a, predicate, root)
    p2_closure = get_closure(graph, pheno_b, predicate, root)
    overlap = p1_closure.intersection(p2_closure)
    max_ic = max([ic_map[parent]for parent in overlap])
    mica = ''
    for pheno in overlap:
        if ic_map[pheno] == max_ic:
            mica = pheno
    return mica


def pairwise_jaccard(pheno_a: str, pheno_b: str, graph: Graph, root:str) -> float:
    predicate = RDFS['subClassOf']
    return math_utils.jaccard(
        get_closure(graph, pheno_a, predicate, root),
        get_closure(graph, pheno_b, predicate, root)
    )


def profile_jaccard(
        profile_a: Sequence[str],
        profile_b: Sequence[str],
        graph: Graph,
        root: str) -> float:
    predicate = RDFS['subClassOf']
    pheno_a_set = {get_closure(graph, pheno, predicate, root) for pheno in profile_a}
    pheno_b_set = {get_closure(graph, pheno, predicate, root) for pheno in profile_b}
    return math_utils.jaccard(pheno_a_set, pheno_b_set)


def hp_iri2curie(iri:str) -> str:
    return iri.replace("http://purl.obolibrary.org/obo/HP_", "HP:")


def hp_curie2iri(curie:str) -> URIRef:
    return URIRef(curie.replace("HP:", "http://purl.obolibrary.org/obo/HP_"))


