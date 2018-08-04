from typing import Sequence, Dict, List, Union, Optional
from enum import Enum
from rdflib.graph import Graph
from rdflib.namespace import RDFS
from phenom.utils import owl_utils, math_utils

# Union types
Num = Union[int, float]


class SimMetric(Enum):
    PHENODIGM = 'phenodigm'
    IC = 'ic'


class Phenodigm():
    """
    Phenodigm algorithm:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3649640/pdf/bat025.pdf

    There are two common variations of the pairwise similarity calculations
    1. geometric mean of the jaccard index and ic of the MICA
    2. ic of the MICA

    The first is the metric used in the published algorithm, the second
    is used in the owltools OWLTools-Sim package
    """

    def __init__(
            self,
            graph: Graph,
            root: str,
            ic_map: Dict[str, float],
            is_same_species: Optional[bool] = True,
            similarity_type: Union[SimMetric, str, None]= SimMetric.PHENODIGM):

        self.graph = graph
        self.root = root
        self.ic_map = ic_map
        self.is_same_species =  is_same_species
        if isinstance(similarity_type, SimMetric):
            self.similarity_type = similarity_type
        else:
            self.similarity_type = SimMetric(similarity_type.lower())

    def phenodigm_compare(
            self,
            profile_a: Sequence[str],
            profile_b: Sequence[str]) -> float:

        query_matrix = self.get_score_matrix(profile_a, profile_b)
        if self.is_same_species:
            optimal_matrix = self.get_optimal_matrix(profile_a)
        else:
            raise NotImplementedError
        return self.compute_phenodigm_score(query_matrix, optimal_matrix)

    def symmetric_phenodigm(
            self,
            profile_a: Sequence[str],
            profile_b: Sequence[str]) -> float:

        a2b_matrix = self.get_score_matrix(list(profile_a), list(profile_b))
        b2a_matrix = math_utils.flip_matrix(a2b_matrix)
        optimal_a_matrix = self.get_optimal_matrix(profile_a)
        optimal_b_matrix = self.get_optimal_matrix(profile_b)

        return math_utils.mean(
            [self.compute_phenodigm_score(a2b_matrix, optimal_a_matrix),
             self.compute_phenodigm_score(b2a_matrix, optimal_b_matrix)])

    @staticmethod
    def compute_phenodigm_score(
            query_matrix: List[List[float]],
            optimal_matrix: List[List[float]]) -> float:
        return 100 * math_utils.mean(
            [math_utils.max_percentage_score(query_matrix, optimal_matrix),
             math_utils.avg_percentage_score(query_matrix, optimal_matrix)])

    def get_mica_ic(self, pheno_a: str, pheno_b: str) -> Num:
        predicate = RDFS['subClassOf']
        p1_closure = owl_utils.get_closure(self.graph, pheno_a, predicate, self.root)
        p2_closure = owl_utils.get_closure(self.graph, pheno_b, predicate, self.root)
        return max([self.ic_map[parent]
                        for parent in p1_closure.intersection(p2_closure)
                    ])

    def jac_ic_geomean(
            self,
            pheno_a: str,
            pheno_b: str) -> float:
        """
        may make more sense in owl_utils
        """
        jaccard_sim = owl_utils.pairwise_jaccard(
            pheno_a, pheno_b, self.graph, self.root)
        mica = self.get_mica_ic(pheno_a, pheno_b)
        return math_utils.geometric_mean([jaccard_sim, mica])

    def get_score_matrix(
            self,
            profile_a: List[str],
            profile_b: List[str]) -> List[List[float]]:

        score_matrix = [[]]

        if self.similarity_type == SimMetric.PHENODIGM:
            sim_fn = self.jac_ic_geomean
        elif self.similarity_type == SimMetric.IC:
            sim_fn = self.get_mica_ic
        else:
            raise NotImplementedError

        for index, pheno_a in enumerate(profile_a):
            if index == len(score_matrix):
                score_matrix.append([])
            for pheno_b in profile_b:
                score_matrix[index].append(
                    sim_fn(pheno_a, pheno_b)
                )
        return score_matrix

    def get_optimal_matrix(self, profile: Sequence[str]) -> List[List[float]]:
        """
        Only implemented for same species comparisons
        """
        score_matrix = []
        if self.is_same_species:
            for pheno in profile:
                if self.similarity_type == SimMetric.PHENODIGM:
                    score_matrix.append(
                        [math_utils.geometric_mean([1, self.ic_map[pheno]])])
                elif self.similarity_type == SimMetric.IC:
                    score_matrix.append([self.ic_map[pheno]])
                else:
                    raise NotImplementedError
        else:
            raise NotImplementedError
        return score_matrix
