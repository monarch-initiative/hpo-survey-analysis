from typing import Sequence, Dict, List, Union, Optional
from enum import Enum
from rdflib import Graph, URIRef, RDFS
from phenom.similarity import metric
from phenom.math import matrix, math
import math as std_math
from functools import reduce


# Union types
Num = Union[int, float]


class PairwiseSim(Enum):
    GEOMETRIC = 'geometric'
    IC        = 'ic'
    JACCARD   = 'jaccard'


class MatrixMetric(Enum):
    MAX = 'max'  # Max
    AVG = 'avg'  # Average
    BMA = 'bma'  # Best Match Average


class SemanticSim():

    def __init__(
            self,
            graph: Graph,
            root: str,
            ic_map: Dict[str, float]):
        self.graph = graph
        self.root = root
        self.ic_map = ic_map

    def sim_gic(
            self,
            profile_a: Sequence[str],
            profile_b: Sequence[str],
            predicate: Optional[URIRef] = RDFS['subClassOf']) -> float:
        """
        Groupwise resnik similarity:
        Summed information content of common ancestors divided by summed
        information content of all ancestors in profile a and profile b
        https://bmcbioinformatics.biomedcentral.com/track/
        pdf/10.1186/1471-2105-9-S5-S4
        """
        a_closure = metric.get_profile_closure(
            profile_a, self.graph, self.root, predicate)
        b_closure = metric.get_profile_closure(
            profile_b, self.graph, self.root, predicate)

        numerator = reduce(
            lambda x, y: x + y,
            [self.ic_map[pheno] for pheno in a_closure.intersection(b_closure)]
        )
        denominator = reduce(
            lambda x, y: x + y,
            [self.ic_map[pheno] for pheno in a_closure.union(b_closure)]
        )

        return numerator/denominator

    def cosine_sim(
            self,
            profile_a: Sequence[str],
            profile_b: Sequence[str],
            predicate: Optional[URIRef] = RDFS['subClassOf']) -> float:
        """
        Implemented as ochai coefficient, or
        len( A union B) / sqrt(len(A)*len(B))

        :param profile_a:
        :param profile_b:
        :return:
        """
        pheno_a_set = metric.get_profile_closure(
            profile_a, self.graph, self.root, predicate)
        pheno_b_set = metric.get_profile_closure(
            profile_b, self.graph, self.root, predicate)
        numerator = len(pheno_a_set.intersection(pheno_b_set))
        denominator = std_math.sqrt(len(pheno_a_set) * len(pheno_b_set))
        return numerator / denominator

    def sim_gicosine(
            self,
            profile_a: Sequence[str],
            profile_b: Sequence[str],
            predicate: Optional[URIRef] = RDFS['subClassOf']) -> float:
        """
        simGICosine
        Cosine similarity where instead of vectors of 0 and 1 for
        present/absent, we use vectors of information content values
        for present classes, and 0 for absent classes
        """
        a_closure = metric.get_profile_closure(
            profile_a, self.graph, self.root, predicate)
        b_closure = metric.get_profile_closure(
            profile_b, self.graph, self.root, predicate)
        numerator = reduce(
            lambda x, y: x + y,
            [std_math.pow(self.ic_map[pheno], 2) for pheno in a_closure.intersection(b_closure)]
        )
        denominator = std_math.sqrt(reduce(
            lambda x, y: x + y,
            [std_math.pow(self.ic_map[pheno], 2) for pheno in a_closure]
        )) * std_math.sqrt(reduce(
            lambda x, y: x + y,
            [std_math.pow(self.ic_map[pheno], 2) for pheno in b_closure]
        ))
        return numerator / denominator

    def jaccard_sim(
            self,
            profile_a: Sequence[str],
            profile_b: Sequence[str],
            predicate: Optional[URIRef] = RDFS['subClassOf']) -> float:
        """
        Groupwise jaccard similarty
        """
        pheno_a_set = metric.get_profile_closure(
            profile_a, self.graph, self.root, predicate)
        pheno_b_set = metric.get_profile_closure(
            profile_b, self.graph, self.root, predicate)

        return metric.jaccard(pheno_a_set, pheno_b_set)

    def resnik_sim(
            self,
            profile_a: Sequence[str],
            profile_b: Sequence[str],
            matrix_metric: Union[MatrixMetric, str, None] = MatrixMetric.BMA,
            is_symmetric: Optional[bool]=False,
            is_normalized: Optional[bool]=False) -> float:
        """
        Resnik similarity

        Implementations:
        BMA (Best Match Average) - Average of the best match per row in matrix
        Max - Max of matrix
        Average - Average IC of entire matrix
        :param profile_a: Sequence of phenotypes
        :param profile_b: Sequence of phenotypes
        :param matrix_metric: BMA, max, avg Default: BMA
        :param is_symmetric: avg(Pa vs Pb, Pb vs Pa) Default: False
        :param is_normalized: Normalize by dividing by the resnik
                              of the optimal matrix Default: False
        :return: resnik score, a float between 0-MaxIC in cache,
                 if normalized a float between 0-1
        """
        if not isinstance(matrix_metric, MatrixMetric):
            matrix_metric = MatrixMetric(matrix_metric.lower())

        resnik_score = 0
        if is_symmetric:
            resnik_score = math.mean(
                [self._compute_resnik_score(profile_a, profile_b,
                                            matrix_metric, is_normalized),
                 self._compute_resnik_score(profile_b, profile_a,
                                            matrix_metric, is_normalized)])
        else:
            resnik_score = self._compute_resnik_score(
                profile_a, profile_b,matrix_metric, is_normalized)

        return resnik_score

    def _compute_resnik_score(
            self,
            profile_a: Sequence[str],
            profile_b: Sequence[str],
            matrix_metric: MatrixMetric,
            is_normalized: bool = False)-> float:

        similarity_type = PairwiseSim.IC

        query_matrix = self._get_score_matrix(
            profile_a, profile_b, similarity_type)

        if is_normalized:
            optimal_matrix = self._get_optimal_matrix(
                profile_a, similarity_type=similarity_type)

        resnik_score = 0

        if matrix_metric == MatrixMetric.BMA:
            if is_normalized:
                resnik_score = matrix.bma_percentage_score(
                    query_matrix, optimal_matrix)
            else:
                resnik_score = matrix.bma_score(query_matrix)
        elif matrix_metric == MatrixMetric.MAX:
            if is_normalized:
                resnik_score = matrix.max_percentage_score(
                    query_matrix, optimal_matrix)
            else:
                resnik_score = matrix.max_score(query_matrix)
        elif matrix_metric == MatrixMetric.AVG:
            if is_normalized:
                resnik_score = matrix.avg_percentage_score(
                    query_matrix, optimal_matrix)
            else:
                resnik_score = matrix.avg_score(query_matrix)

        return resnik_score

    def phenodigm_compare(
            self,
            profile_a: Sequence[str],
            profile_b: Sequence[str],
            is_symmetric: Optional[bool]=False,
            is_same_species: Optional[bool]=True,
            similarity_type: Union[PairwiseSim, str, None]= PairwiseSim.GEOMETRIC
    ) -> float:
        """
        Phenodigm algorithm:
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3649640/pdf/bat025.pdf

        There are two common variations of the pairwise similarity calculations
        1. geometric mean of the jaccard index and ic of the MICA
        2. ic of the MICA

        The first is the metric used in the published algorithm, the second
        is used in the owltools OWLTools-Sim package
        """
        if not isinstance(similarity_type, PairwiseSim):
            similarity_type = PairwiseSim(similarity_type.lower())

        query_matrix = self._get_score_matrix(
            profile_a, profile_b, similarity_type)
        if is_same_species:
            optimal_matrix = self._get_optimal_matrix(
                profile_a, is_same_species, similarity_type)
        else:
            raise NotImplementedError

        if is_symmetric:
            b2a_matrix = matrix.flip_matrix(query_matrix)
            optimal_b_matrix = self._get_optimal_matrix(
                profile_b, is_same_species, similarity_type)
            score = math.mean(
                [self.compute_phenodigm_score(query_matrix, optimal_matrix),
                 self.compute_phenodigm_score(b2a_matrix, optimal_b_matrix)])
        else:
            score = self.compute_phenodigm_score(query_matrix, optimal_matrix)

        return score

    @staticmethod
    def compute_phenodigm_score(
            query_matrix: List[List[float]],
            optimal_matrix: List[List[float]]) -> float:
        return 100 * math.mean(
            [matrix.max_percentage_score(query_matrix, optimal_matrix),
             matrix.bma_percentage_score(query_matrix, optimal_matrix)])

    def _get_score_matrix(
            self,
            profile_a: Sequence[str],
            profile_b: Sequence[str],
            similarity_type:Union[PairwiseSim, None] = PairwiseSim.IC
    ) -> List[List[float]]:

        score_matrix = [[]]

        if similarity_type == PairwiseSim.GEOMETRIC:
            sim_fn = metric.jac_ic_geomean
        elif similarity_type == PairwiseSim.IC:
            sim_fn = metric.get_mica_ic
        else:
            raise NotImplementedError

        for index, pheno_a in enumerate(profile_a):
            if index == len(score_matrix):
                score_matrix.append([])
            for pheno_b in profile_b:
                score_matrix[index].append(
                    sim_fn(pheno_a, pheno_b, self.graph, self.ic_map, self.root)
                )
        return score_matrix

    def _get_optimal_matrix(
            self,
            profile: Sequence[str],
            is_same_species: Optional[bool]=True,
            similarity_type: Union[PairwiseSim, None]= PairwiseSim.IC
    ) -> List[List[float]]:
        """
        Only implemented for same species comparisons
        """
        score_matrix = []
        if is_same_species:
            for pheno in profile:
                if similarity_type == PairwiseSim.GEOMETRIC:
                    score_matrix.append(
                        [math.geometric_mean([1, self.ic_map[pheno]])])
                elif similarity_type == PairwiseSim.IC:
                    score_matrix.append([self.ic_map[pheno]])
                else:
                    raise NotImplementedError
        else:
            raise NotImplementedError
        return score_matrix
