from typing import Dict, Set, FrozenSet, Optional
from phenom.utils.owl_utils import get_closure
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
    omission_rate = .2
    imprecision_rate = .1
    noise_rate = .3

    phenotypes = list(pheno_profile)
    profile_size = len(phenotypes)

    # Remove 20 percent of phenotypes
    count_to_remove = round(profile_size * omission_rate)
    phenotypes = random.sample(phenotypes, profile_size - count_to_remove)

    # mutate 10 percent to closest parent
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
