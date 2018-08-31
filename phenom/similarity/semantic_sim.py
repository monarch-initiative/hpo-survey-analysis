from typing import Iterable, Dict, List, Union, Optional
from enum import Enum
from rdflib import Graph, URIRef, RDFS
from phenom.similarity import metric
from phenom.math import matrix, math_utils
import math
from functools import reduce
import numpy as np


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
            profile_a: Iterable[str],
            profile_b: Iterable[str],
            predicate: Optional[URIRef] = RDFS['subClassOf']) -> float:
        """
        Groupwise resnik similarity:
        Summed information content of common ancestors divided by summed
        information content of all ancestors in profile a and profile b
        https://bmcbioinformatics.biomedcentral.com/track/
        pdf/10.1186/1471-2105-9-S5-S4
        """
        # Filter out negative phenotypes
        profile_a = {pheno for pheno in profile_a if not pheno.startswith("-")}
        profile_b = {pheno for pheno in profile_b if not pheno.startswith("-")}

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
            profile_a: Iterable[str],
            profile_b: Iterable[str],
            ic_weighted: Optional[bool] = False,
            negative_weight: Optional[Num] = 1,
            predicate: Optional[URIRef] = RDFS['subClassOf']) -> float:
        """
        Cosine similarity
        Profiles are treated as vectors of numbers between 0-1:
        1: Phenotype present
        0: Absent (no information)
        1 * negative weight: Negated phenotypes

        if ic_weighted is true the attributes become vectors
        of information content scores

        Inferred phenotypes are computed as parent classes for positive phenotypes
        and child classes for negative phenotypes.  Typically we do not want to
        weight negative phenotypes as high as positive phenotypes.  A weight between
        .01-.1 may be desirable
        """
        def score(term):
            if ic_weighted:
                attribute = self.ic_map[term]
            else:
                attribute = 1
            return attribute

        positive_a_profile = {item for item in profile_a if not item.startswith('-')}
        negative_a_profile = {item[1:] for item in profile_a if item.startswith('-')}

        positive_b_profile = {item for item in profile_b if not item.startswith('-')}
        negative_b_profile = {item[1:] for item in profile_b if item.startswith('-')}

        pos_a_closure = metric.get_profile_closure(
            positive_a_profile, self.graph, self.root, predicate)
        pos_b_closure = metric.get_profile_closure(
            positive_b_profile, self.graph, self.root, predicate)

        neg_a_closure = {"-{}".format(item)
                              for item in metric.get_profile_closure(
            negative_a_profile, self.graph, self.root, predicate, negative=True)
        }

        neg_b_closure = {"-{}".format(item)
                              for item in metric.get_profile_closure(
            negative_b_profile, self.graph, self.root, predicate, negative=True)
        }

        pos_intersect_dot_product = reduce (
            lambda x, y: x + y,
            [math.pow(score(item), 2)
             for item in pos_a_closure.intersection(pos_b_closure)],
            0
        )

        neg_intersect_dot_product = reduce(
            lambda x, y: x + y,
            [math.pow(score(item) * negative_weight, 2)
             for item in neg_a_closure.intersection(neg_b_closure)],
            0
        )

        a_square_dot_product = math.sqrt(
            reduce(
                lambda x, y: x + y,
                [math.pow(score(item), 2) for item in pos_a_closure],
                0
            ) +
            reduce(
                lambda x, y: x + y,
                [math.pow(score(item) * negative_weight, 2)
                 for item in neg_a_closure],
                0
            )
        )

        b_square_dot_product = math.sqrt(
            reduce(
                lambda x, y: x + y,
                [math.pow(score(item), 2) for item in pos_b_closure],
                0
            ) +
            reduce(
                lambda x, y: x + y,
                [math.pow(score(item) * negative_weight, 2)
                 for item in neg_b_closure],
                0
            )
        )

        numerator = pos_intersect_dot_product + neg_intersect_dot_product
        denominator = a_square_dot_product * b_square_dot_product

        return numerator / denominator

    def jaccard_sim(
            self,
            profile_a: Iterable[str],
            profile_b: Iterable[str],
            predicate: Optional[URIRef] = RDFS['subClassOf']) -> float:
        """
        Groupwise jaccard similarty
        Negative phenotypes must be prefixed with a '-'
        """
        # Filter out negative phenotypes
        profile_a = {pheno for pheno in profile_a if not pheno.startswith("-")}
        profile_b = {pheno for pheno in profile_b if not pheno.startswith("-")}

        pheno_a_set = metric.get_profile_closure(
            profile_a, self.graph, self.root, predicate)
        pheno_b_set = metric.get_profile_closure(
            profile_b, self.graph, self.root, predicate)

        return metric.jaccard(pheno_a_set, pheno_b_set)

    def resnik_sim(
            self,
            profile_a: Iterable[str],
            profile_b: Iterable[str],
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
        # Filter out negative phenotypes
        profile_a = {pheno for pheno in profile_a if not pheno.startswith("-")}
        profile_b = {pheno for pheno in profile_b if not pheno.startswith("-")}

        if not isinstance(matrix_metric, MatrixMetric):
            matrix_metric = MatrixMetric(matrix_metric.lower())

        sim_measure = PairwiseSim.IC

        query_matrix = self._get_score_matrix(
            profile_a, profile_b, sim_measure)

        if is_normalized:
            optimal_matrix = self._get_optimal_matrix(
                profile_a, sim_measure=sim_measure)
        else:
            optimal_matrix = None

        resnik_score = 0
        if is_symmetric:
            b2a_matrix = matrix.flip_matrix(query_matrix)
            optimal_b_matrix = self._get_optimal_matrix(
                profile_b, sim_measure=sim_measure)
            resnik_score = math_utils.mean(
                [self._compute_resnik_score(
                    query_matrix, optimal_matrix, matrix_metric),
                 self._compute_resnik_score(
                     b2a_matrix, optimal_b_matrix, matrix_metric)])
        else:
            resnik_score = self._compute_resnik_score(
                query_matrix, optimal_matrix, matrix_metric)

        return resnik_score

    def _compute_resnik_score(
            self,
            query_matrix: List[List[float]],
            optimal_matrix: Optional[List[List[float]]] = None,
            matrix_metric: Optional[MatrixMetric] = MatrixMetric.BMA )-> float:

        is_normalized = True if optimal_matrix else False

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
            profile_a: Iterable[str],
            profile_b: Iterable[str],
            is_symmetric: Optional[bool]=False,
            is_same_species: Optional[bool]=True,
            sim_measure: Union[PairwiseSim, str, None]= PairwiseSim.GEOMETRIC
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
        # Filter out negative phenotypes
        profile_a = {pheno for pheno in profile_a if not pheno.startswith("-")}
        profile_b = {pheno for pheno in profile_b if not pheno.startswith("-")}

        if not isinstance(sim_measure, PairwiseSim):
            sim_measure = PairwiseSim(sim_measure.lower())

        query_matrix = self._get_score_matrix(profile_a, profile_b, sim_measure)
        if is_same_species:
            optimal_matrix = self._get_optimal_matrix(
                profile_a, is_same_species, sim_measure)
        else:
            raise NotImplementedError

        if is_symmetric:
            b2a_matrix = matrix.flip_matrix(query_matrix)
            optimal_b_matrix = self._get_optimal_matrix(
                profile_b, is_same_species, sim_measure)
            score = math_utils.mean(
                [self.compute_phenodigm_score(query_matrix, optimal_matrix),
                 self.compute_phenodigm_score(b2a_matrix, optimal_b_matrix)])
        else:
            score = self.compute_phenodigm_score(query_matrix, optimal_matrix)

        return score

    @staticmethod
    def compute_phenodigm_score(
            query_matrix: List[List[float]],
            optimal_matrix: List[List[float]]) -> float:
        return 100 * math_utils.mean(
            [matrix.max_percentage_score(query_matrix, optimal_matrix),
             matrix.bma_percentage_score(query_matrix, optimal_matrix)])

    def _get_score_matrix(
            self,
            profile_a: Iterable[str],
            profile_b: Iterable[str],
            sim_measure: Union[PairwiseSim, None] = PairwiseSim.IC
    ) -> List[List[float]]:

        score_matrix = [[]]

        if sim_measure == PairwiseSim.GEOMETRIC:
            sim_fn = metric.jac_ic_geomean
        elif sim_measure == PairwiseSim.IC:
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
            profile: Iterable[str],
            is_same_species: Optional[bool]=True,
            sim_measure: Union[PairwiseSim, None]= PairwiseSim.IC
    ) -> List[List[float]]:
        """
        Only implemented for same species comparisons
        """
        score_matrix = []
        if is_same_species:
            for pheno in profile:
                if sim_measure == PairwiseSim.GEOMETRIC:
                    score_matrix.append(
                        [math_utils.geometric_mean([1, self.ic_map[pheno]])])
                elif sim_measure == PairwiseSim.IC:
                    score_matrix.append([self.ic_map[pheno]])
                else:
                    raise NotImplementedError
        else:
            raise NotImplementedError
        return score_matrix
