from typing import Dict, Set, FrozenSet, Optional, List
from phenom.utils.owl_utils import get_closure
from phenom.math.math_utils import binomial_coeff
from rdflib import Graph, RDFS
import random
import logging


def simulate_from_derived(
        pheno_profile: Set[str],
        pheno_subset: Set[str],
        graph: Graph,
        root: str,
        ic_values: Dict[str, float],
        filter_out: Set[str],
        ref_disease: Optional[str]=None) -> FrozenSet[str]:
    """
    Add imprecision and noise to profile
    20% omit phenotype - omissions
    10% use closest parent - imprecision
    30% add random phenotype - noise, min 1
    :return: FrozenSet[str] - set of phenotype curies
    """
    omission_rate = .4  # .4 for gold, .2 for derived
    imprecision_rate = .3  # .3 for gold, .1 for derived
    noise_rate = .3

    phenotypes = list(pheno_profile)
    profile_size = len(phenotypes)

    # Remove x percent of phenotypes
    count_to_remove = round(profile_size * omission_rate)
    phenotypes = random.sample(phenotypes, profile_size - count_to_remove)

    # mutate x percent to closest parent
    count_to_mutate = round(profile_size * imprecision_rate)
    random.shuffle(phenotypes)
    counter = 0
    for idx, pheno in enumerate(phenotypes):
        if counter == count_to_mutate:
            break
        parents = get_closure(graph, pheno, RDFS['subClassOf'], root, False)
        lay_overlap = parents.intersection(pheno_subset).difference(pheno_profile, phenotypes)
        if len(list(lay_overlap)) == 0:
            continue
        max_ic = max([ic_values[parent] for parent in lay_overlap])
        mica = ''
        for phen in lay_overlap:
            if ic_values[phen] == max_ic:
                mica = phen

        phenotypes[idx] = mica
        counter += 1

    if counter != count_to_mutate:
        logging.info("Could not mutate profile derived from {}".format(ref_disease))

    # add random phenotype(s)
    # Filter out phenotypes from filter_out set
    phenos_to_select = pheno_subset.difference(filter_out, phenotypes, pheno_profile)

    if len(list(phenos_to_select)) == 0:
        logging.warning("No phenotypes to select for "
                        "profile derived from {}".format(ref_disease))
    comissions = round(profile_size * noise_rate)
    noise_count = 1 if comissions == 0 else comissions

    for i in range(noise_count):
        random_pheno = random.choice(list(phenos_to_select))
        phenotypes.append(random_pheno)
        phenos_to_select.remove(random_pheno)

    return frozenset(phenotypes)


def average_ties(previous_rank: int, tie_count: int) -> int:
    deranked_summed = \
        binomial_coeff(previous_rank + (tie_count)) - \
        binomial_coeff(previous_rank)
    return round(deranked_summed / tie_count)


def rerank_by_average(ranks: List[int]) -> List[int]:
    """
    Given a list of ranked results, averages ties and
    reranks following
    :param matches:
    :return: List of ranks
    """
    # list to store adjusted ranks
    avg_ranks = []
    # Resolve ties, take the average rank
    last_rank = 1
    last_avg_rank = 0
    tie_count = 1
    is_first = True
    for rank in ranks[1:]:
        if rank > last_rank:
            if tie_count == 1:
                last_avg_rank += 1
                if is_first:
                    avg_ranks.append(1)
                    is_first = False
                else:
                    avg_ranks.append(last_avg_rank)
            else:
                avg_rank = average_ties(last_avg_rank, tie_count)
                tied_ranks = [avg_rank for n in range(tie_count)]
                avg_ranks.extend(tied_ranks)
                # reset tie count
                last_avg_rank = avg_rank
                tie_count = 1
        else:
            tie_count += 1
            is_first = False

        last_rank = rank

    if tie_count > 0:
        avg_rank = average_ties(last_avg_rank, tie_count)
        tied_ranks = [avg_rank for n in range(tie_count)]
        avg_ranks.extend(tied_ranks)
    else:
        last_avg_rank += 1
        avg_ranks.append(last_avg_rank)

    return avg_ranks
