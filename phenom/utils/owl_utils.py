from typing import Set, List, Optional, Dict, Sequence
from phenom.decorators import memoized
from rdflib import URIRef, BNode, Literal, Graph, RDFS
from prefixcommons import contract_uri, expand_uri


def get_closure(
        graph: Graph,
        node: str,
        edge: Optional[URIRef]=None,
        root: Optional[str]=None,
        reflexive: Optional[bool] = True,
        negative: Optional[bool] = False) -> Set[str]:

    prefix, reference = node.split(":")
    reference = int(reference)
    closure = get_closure_memoized(
        graph, reference, prefix, edge, root, reflexive, negative)
    return {"{}:{}".format(prefix, str(node).zfill(7)) for node in closure}


def _get_closure(
        graph: Graph,
        node: str,
        edge: Optional[URIRef]=None,
        root: Optional[str]=None,
        reflexive: Optional[bool] = True,
        negative: Optional[bool] = False) -> Set[str]:
    nodes = set()
    if negative:
        nodes = get_descendants(graph, node, edge)
    else:
        nodes = get_ancestors(graph, node, edge, root, reflexive)
    return nodes

@memoized
def get_closure_memoized(
        graph: Graph,
        node_id: int,
        curie_prefix: str,
        edge: Optional[URIRef]=None,
        root: Optional[str]=None,
        reflexive: Optional[bool] = True,
        negative: Optional[bool] = False) -> List[int]:

    node = "{}:{}".format(curie_prefix, str(node_id).zfill(7))
    closure = _get_closure(graph, node, edge, root, reflexive, negative)
    return [int(id.replace("{}:".format(curie_prefix), "")) for id in closure]


def get_mica_id(
        pheno_a: str,
        pheno_b: str,
        graph: Graph,
        ic_map:Dict[str, float],
        root: str) -> str:
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


def get_descendants(
        graph: Graph,
        node: str,
        edge: Optional[URIRef]=RDFS['subClassOf'],
        reflexive: Optional[bool] = True) -> Set[str]:

    nodes = set()
    node = URIRef(expand_uri(node, strict=True))
    for sub in graph.transitive_subjects(edge, node):
        if not reflexive and node == sub:
            continue
        if isinstance(sub, Literal):
            continue
        nodes.add(contract_uri(str(sub), strict=True)[0])
    return nodes


def get_ancestors(
        graph: Graph,
        node: str,
        edge: Optional[URIRef] = RDFS['subClassOf'],
        root: Optional[str] = None,
        reflexive: Optional[bool] = True) -> Set[str]:
    nodes = set()
    root_seen = {}
    if node is not None:
        node = URIRef(expand_uri(node, strict=True))
    if root is not None:
        root = URIRef(expand_uri(root, strict=True))
        root_seen = {root: 1}
    for obj in graph.transitive_objects(node, edge, root_seen):
        if isinstance(obj, Literal) or isinstance(obj, BNode):
            continue
        if not reflexive and node == obj:
            continue
        nodes.add(contract_uri(str(obj), strict=True)[0])

    # Add root to graph
    if root is not None:
        nodes.add(contract_uri(str(root), strict=True)[0])

    return nodes
