from typing import Dict, Set, FrozenSet, Optional, List
from phenom.utils.owl_utils import get_closure
from phenom.math.math_utils import binomial_coeff
from phenom.monarch import owlsim_classify
from phenom.model.synthetic import SyntheticProfile
from rdflib import Graph, RDFS
from multiprocessing import Queue
import numpy
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
    omission_rate = .2  # .4 for gold, .2 for derived
    imprecision_rate = .1  # .3 for gold, .1 for derived
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


def create_confusion_matrix_per_threshold(
        synthetic_profiles: List[SyntheticProfile],
        owlsim_url: str,
        num_labels: int,
        thresholds: List,
        threshold_type: str):

    confusion_by_rank: Dict[int, List[int]] = {}

    counter = 0
    total = len(synthetic_profiles)

    for threshold in thresholds:
        confusion_by_rank[threshold] = [0, 0, 0, 0]

    for synth_profile in synthetic_profiles:

        if counter % 1000 == 0:
            logging.info("processed {} patients out of {}".format(counter, total))
        counter += 1

        sim_resp = owlsim_classify(synth_profile.phenotypes, num_labels, owlsim_url)

        if 'matches' not in sim_resp:
            raise ValueError("Could not find match for {}".format(synth_profile.disease))

        # could dynamically generated num_classes here
        # num_classes = len(sim_resp['matches'])

        # store the disease rank or score
        disease_score = 0
        for disease_index, match in enumerate(sim_resp['matches']):
            if match['matchId'] == synth_profile.disease:
                if threshold_type == 'rank':
                    disease_score = disease_index
                elif threshold_type == 'probability':
                    disease_score = float(match['rawScore'])

        # Create a list to track the number of positive hits per threshold
        # This list should have the same number of indices as thresholds,
        # conceptually they are two columns in a table

        # For rank, we count the number of diseases for each rank, and compute
        # the number of diseases per rank threshold by using sum,
        # eg number of diseases <= rank 4 is sum(positives[0:5])

        # For probability, we count the number of labels <= to each probability
        # therefore to get the total count at an index, get the val of positives[index]
        positives = []

        if threshold_type == 'rank':
            positives = [0 for r in range(0, len(thresholds) + 1)]
            owlsim_ranks = [match['rank'] for match in sim_resp['matches']]
            avg_ranks = rerank_by_average(owlsim_ranks)
            disease_score = avg_ranks[disease_score]

            for rnk in avg_ranks:
                positives[rnk] += 1
            positives.pop(0)  # ranks start at 1, so make 1st rank the 0th index

        elif threshold_type == 'probability':
            scores = numpy.array([float(match['rawScore']) for match in sim_resp['matches']])
            positives = [0 for r in range(0, len(thresholds))]
            for index, prob in enumerate(thresholds):
                score_count = 0
                if index == 0:
                    positives[index] = 0
                else:
                    positives[index] = positives[index-1]
                for score in scores:
                    if score >= prob:
                        score_count += 1
                    else:
                        break
                positives[index] += score_count
                scores = scores[score_count:]

        else:
            raise ValueError("{} not valid threshold".format(threshold_type))

        for index, threshold in enumerate(thresholds):
            true_pos, false_pos, false_neg, true_neg = confusion_by_rank[threshold]

            is_true_pos = False
            if threshold_type == 'rank':
                is_true_pos = disease_score <= threshold
                positive_diseases = sum(positives[0:threshold])
            elif threshold_type == 'probability':
                is_true_pos = disease_score >= threshold
                positive_diseases = positives[index]
            else:
                raise ValueError

            if is_true_pos:
                true_pos += 1
                true_neg += num_labels - positive_diseases
                false_pos += positive_diseases - 1
            else:
                false_neg += 1
                false_pos += positive_diseases
                true_neg += num_labels - positive_diseases - 1

            confusion_by_rank[threshold] = [true_pos, false_pos, false_neg, true_neg]

    return confusion_by_rank


def process_confusion_matrix_per_threshold(
        synthetic_profiles: List[SyntheticProfile],
        owlsim_url: str,
        num_labels: int,
        thresholds: List,
        threshold_type: str,
        queue: Queue) -> None:
    """
    Run create_confusion_matrix_per_threshold and put results in queue object
    """
    queue.put(create_confusion_matrix_per_threshold(
        synthetic_profiles,
        owlsim_url,
        num_labels,
        thresholds,
        threshold_type
    ))
