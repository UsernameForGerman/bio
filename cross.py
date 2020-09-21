from typing import Union, List
import numpy as np
from random import choices

from .models import Animal, Bull, Cow, Genotype, Phenotype

def get_genotype_indexies(r: List[float]) -> List[int]:
    """
    Get indexies for choising haplotype in genotype
    """
    return_list = list()
    for i in range(1, len(r) + 1):
        if i == 1:
            value = choices((0, 1), weights=(0.5, 0.5))[0]
        else:
            value = choices((return_list[i - 1], 1 - return_list[i - 1]), weights=(1 - r[i - 1], r[i - 1]))[0]

        return_list.append(value)

    return return_list


def get_gamete(genotype: List[List[int]], indexies: List[int]) -> List[int]:
    """
    Return gamete from genotype input
    """
    return [
        genotype[i][indexies[i]] for i in range(len(genotype[0]))
    ]

def get_reproduce(genotypes: List[List[List[int]]], possibilities: List[float], population_of_progeny: int) -> List[List[List[int]]]:
    """
    Return children(which is amount=population_of_progeny)' genotype from parents' genotypes
    """
    haplotypes = [[[]]]
    for i in range(len(genotypes[0])):
        for j in range(len(genotypes)):
            for k in range(population_of_progeny):
                haplotypes[i][j][k] = get_gamete(genotypes[j], get_genotype_indexies(possibilities))
    return haplotypes

def gebv_select(reproduce: List[List[List[int]]], possibilities: List[float]):
    return

def pcv_select():
    return

def ohv_select():
    return

def wgebv_select():
    return

def ie312_select():
    return

def ie534_select():
    return

def ie634_select():
    return

def select(reproduce: List[List[List[int]]], possibilities: List[float], strategy: str='gebv') -> (int, int):
    """

    Parameters
    ----------
    strategy: str, gebv, pcv, ohv, wgebv, ie312, ie534, ie634

    Returns
    -------

    """
    if strategy == 'gebv':
        return gebv_select(reproduce, possibilities)

