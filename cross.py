from typing import Union, List, Iterable, Dict
import numpy as np
from random import choices

from .models import Animal, Bull, Cow, Genotype, Phenotype


#
#
# def get_gamete(genotype: List[List[int]], indexies: List[int]) -> List[int]:
#     """
#     Return gamete from genotype input
#     """
#     return [
#         genotype[i][indexies[i]] for i in range(len(genotype[0]))
#     ]
#
#
# def get_reproduce(genotypes: np.ndarray, possibilities: List[float], population_of_progeny: int):
#     """
#     Return children(which is amount=population_of_progeny)' genotype from parents' genotypes
#     """
#     haplotypes = [[[]]]
#     for i in range(len(genotypes[0])):
#         for j in range(len(genotypes)):
#             for k in range(population_of_progeny):
#                 haplotypes[i][j][k] = get_gamete(genotypes[j], get_genotype_indexies(possibilities))
#     return np.array(haplotypes)


class Feature:
    @classmethod
    def get_feature(cls, approach="gebv"):
        if approach == "gebv":
            return getattr(cls, "_gebv_feature")

    @staticmethod
    def gebv_feature(genotype) -> int:
        return np.sum(genotype[0, :] + genotype[1, :])


class Breed:
    def __init__(self, parent_genotypes: np.ndarray, possibilites: List[float], population_of_progeny: int,
                 maximum_feature=None):
        self.generations: Dict[int, np.ndarray] = {0: parent_genotypes}
        self.possibilities = possibilites
        self.population_of_progeny = population_of_progeny
        self._current_generation = 0
        self.maximum_feature = self.generations[0].shape[2] * 2 if maximum_feature is None else maximum_feature

    def evaluate(self, selection: str="gebv", max_generations=None):
        if not getattr(self, f"_{selection}_selection", None):
            raise NotImplementedError(f"Селекция {selection} не реализована!")

        current_generation_number = 0
        while True:
            if self._is_max_feature_genotype(current_generation_number):
                return current_generation_number
            selected_parents = getattr(self, f"_{selection}_selection")(current_generation_number)
            parents_genotypes = self.generations[current_generation_number][selected_parents, :, :]
            self.generations[current_generation_number + 1] = self.get_reproduce(parents_genotypes)
            if max_generations and current_generation_number == max_generations:
                return current_generation_number
            current_generation_number += 1

    def _gebv_selection(self, current_generation_number) -> Iterable[int]:
        reproduce = self.generations[current_generation_number]
        population = reproduce.shape[0]
        features = {Feature.gebv_feature(reproduce[i, :, :]): i for i in range(population)}
        sorted_features = sorted(features)
        return features[sorted_features[-1]], features[sorted_features[-2]]

    def get_reproduce(self, parent_genotypes: np.ndarray):
        """
        Return children(which is amount=population_of_progeny)' genotype from parents' genotypes
        """
        haplotypes = np.zeros((self.population_of_progeny, 2, parent_genotypes.shape[2]))
        number_of_parents = parent_genotypes.shape[0]
        for k in range(self.population_of_progeny):
            new_genotype = None
            for j in range(number_of_parents):
                if new_genotype:
                    np.append(
                        new_genotype, self.get_gamete(parent_genotypes[j], self.get_genotype_indexies()),
                        axis=0
                    )
                else:
                    new_genotype = self.get_gamete(parent_genotypes[j], self.get_genotype_indexies())
            haplotypes[k, :, :] = new_genotype
        return haplotypes

    def get_genotype_indexies(self) -> List[int]:
        """
        Get indexies for choising haplotype in genotype
        """
        return_list = list()
        for i in range(1, len(self.possibilities) + 1):
            if i == 1:
                value = choices((0, 1), weights=(0.5, 0.5))[0]
            else:
                value = choices(
                    (return_list[i - 1], 1 - return_list[i - 1]),
                    weights=(1 - self.possibilities[i - 1], self.possibilities[i - 1])
                )[0]

            return_list.append(value)

        return return_list

    @staticmethod
    def get_gamete(genotype: np.ndarray, indexies: List[int]) -> List[int]:
        """
        Return gamete from genotype input
        """
        return [
            genotype[indexies[i], i] for i in range(genotype.shape[1])
        ]

    def _is_max_feature_genotype(self, current_generation_number):
        for genotype in self.generations[current_generation_number]:
            if np.sum(genotype) >= self.maximum_feature:
                return True
        return False

#
# def gebv_select(reproduce: List[List[List[int]]], possibilities: List[int], n: int) -> int:
#     return np.sum(possibilities * np.array(reproduce[n])[:, 1]) + np.sum(possibilities * np.array(reproduce[n])[:, 0])
#
#
# def ohv_select(reproduce: List[List[List[int]]], possibilities: List[int], n: int) -> int:
#     # TODO: maximum takes only 2 args
#     return np.sum(2 * np.maximum(
#         possibilities * np.array(reproduce[n])[:, 1],
#         possibilities * np.array(reproduce[n])[:, 0])
#     )
#
#
# def opv_select():
#     return
#
#
# def pcv_select():
#     return
#
# def gb_select():
#     return
#
#
# def wgebv_select(
#         reproduce: List[List[List[int]]],
#         weight_list: List[int],
#         possibilities: List[int], n: int, ) -> int:
#
#     weighted_possib = possibilities / \
#                       (np.sqrt(np.maximum(np.array([1 / len(reproduce)] *
#                     len(reproduce[n])[:, 0]), weight_list)))
#
#
#     return (np.sum(weighted_possib * np.array(reproduce[n])[:, 1]) +
#            np.sum(weighted_possib * np.array(reproduce[n])[:, 0]))
#
#
# def ie312_select():
#     return
#
#
# def ie534_select():
#     return
#
#
# def ie634_select():
#     return
#
#
# def select(reproduce: List[List[List[int]]], possibilities: List[float], strategy: str = 'gebv') -> (int, int):
#     """
#
#     Parameters
#     ----------
#     strategy: str, gebv, pcv, ohv, wgebv, ie312, ie534, ie634
#
#     Returns
#     -------
#
#     """
#     if strategy == 'gebv':
#         return gebv_select(reproduce, possibilities)
