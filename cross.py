from typing import Union, List, Iterable, Dict, Callable, Any
import numpy as np
from random import choices
from dataclasses import dataclass, fields

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


@dataclass
class DefaultVal:
    val: Any


@dataclass
class DefaultValuesMixin:
    def __post_init__(self):
        for field in fields(self):

            # if a field of this data class defines a default value of type
            # `DefaultVal`, then use its value in case the field after
            # initialization has either not changed or is None.
            if isinstance(field.default, DefaultVal):
                field_val = getattr(self, field.name)
                if isinstance(field_val, DefaultVal) or field_val is None:
                    setattr(self, field.name, field.default.val)


@dataclass
class Genotype:
    matrix: np.ndarray


@dataclass
class Generation(DefaultValuesMixin):
    index: int
    genotypes: List[Genotype]
    population: int = 1


@dataclass
class Feature:
    index: int
    number: int


def filter_generation_by_index(generations: List[Generation], generation_index: int = 0):
    generation_list = list(filter(lambda generation: generation.index == generation_index, generations))
    if not generation_list:
        raise AttributeError(f"Поколения №{generation_index} не существует в переданном массиве")
    return generation_list[0]


class FeatureHelper:
    @classmethod
    def get_feature(cls, approach="gebv"):
        if approach == "gebv":
            return getattr(cls, "_gebv_feature")

    @classmethod
    def _gebv_feature(cls, genotype: Genotype) -> int:
        return np.sum(genotype.matrix[0, :] + genotype.matrix[1, :])

    @classmethod
    def calculate_features_for_generations(cls, generations: List[Generation], approach: str = None,
                                           feature_function: Callable[[Genotype], Feature] = None) -> List[Feature]:
        _feature_function = cls.__get_valid_feature_function(approach, feature_function)
        generation_features = []
        for generation in generations:
            np.append(generation_features, cls._calculate_feature_for_generation(generation, _feature_function), axis=0)
        return generation_features

    @classmethod
    def calculate_feature_for_generation(
            cls, generations: List[Generation], generation_number: int = 0,
            approach: str = None, feature_function: Callable[[Genotype], Feature] = None) -> List[Feature]:
        _feature_function = cls.__get_valid_feature_function(approach, feature_function)
        filter_generation_by_index(generations, generation_number)

        return cls._calculate_feature_for_generation(generation_list[0], _feature_function)

    @classmethod
    def _calculate_feature_for_generation(cls, generation: Generation, feature_function: Callable[[Genotype], Feature]
                                          ) -> List[Feature]:
        generation_features = []
        for genotype in generation.genotypes:
            generation_features.append(feature_function(genotype))
        return generation_features

    @classmethod
    def __get_valid_feature_function(cls, approach, feature_function):
        """
        Return selected feature function
        """
        if not approach and not feature_function:
            raise AttributeError("Выберите метод вычисления или укажите явно функцию вычисления")

        _feature_function = feature_function
        if not _feature_function:
            if not getattr(cls, f"_{approach}_selection", None):
                raise NotImplementedError(f"Метод вычисления {approach} не реализован!")
            _feature_function = getattr(cls, f"_{approach}_selection")

        return _feature_function


class GraphHelper:

    @classmethod
    def _feature_individual(cls, genotype: Genotype, feature_function: Callable[[Genotype], Feature]):
        return

    @classmethod
    def feature_individual_for_generation(cls, generation: Generation, feature_function: Callable[[Genotype], Feature]):
        """
        Plot graphic
        """
        features = []
        indexs = []
        i = 0
        for genotype in generation.genotypes:
            features.append(cls._feature_individual(genotype, feature_function))
            indexs.append(i)
            i += 1

        return features, indexs

    @classmethod
    def pcv_map(cls):
        """
        Plot graphic
        """
        return

    @classmethod
    def proportion_of_desirable_alleles(cls, selection: str = "gebv"):
        """
        Plot graphic
        """
        return


class Breed:
    # TODO: Convert all vars to dataclasses
    def __init__(self, parent_genotypes: np.ndarray, possibilites: List[float], population_of_progeny: int,
                 maximum_feature=None):
        if len(parent_genotypes.shape) != 3 or parent_genotypes.shape[:2] != (2, 2) or parent_genotypes.shape[2] <= 0:
            raise AttributeError("Массив генотипов особей задан неверно! Размерность должна быть (2 x 2 x N)")

        # self.generations = {0: parent_genotypes}
        self.generations: List[Generation] = [Generation(
            index=0,
            genotypes=list(Genotype(genotype) for genotype in parent_genotypes),
            population=parent_genotypes.shape[0]
        )]
        self.possibilities = possibilites
        self.population_of_progeny = population_of_progeny
        self._current_generation = 0
        self.maximum_feature = (self.generations[0].genotypes[0].matrix.shape[2] * 2
                                if maximum_feature is None
                                else maximum_feature)

    def evaluate(self, selection: str = "gebv", max_generations: int = None):
        if not getattr(self, f"_{selection}_selection", None):
            raise NotImplementedError(f"Селекция {selection} не реализована!")

        current_generation_number = 0
        while True:
            if self._is_max_feature_genotype(current_generation_number):
                return current_generation_number
            selected_parents = getattr(FeatureHelper, f"_{selection}_selection")(current_generation_number)
            parents_genotypes = self.generations[current_generation_number][selected_parents, :, :]
            self.generations[current_generation_number + 1] = self.get_reproduce(parents_genotypes)
            if max_generations and current_generation_number == max_generations:
                return current_generation_number
            current_generation_number += 1

    def _gebv_selection(self, current_generation_number) -> Iterable[int]:
        reproduce = self.generations[current_generation_number]
        population = reproduce.shape[0]
        features = {FeatureHelper._gebv_feature(reproduce[i, :, :]): i for i in range(population)}
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
        return_list = []
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

    def _is_max_feature_genotype(self, current_generation_number: int):
        generation = filter_generation_by_index(self.generations, current_generation_number)
        for genotype in generation.genotypes:
            if np.sum(genotype.matrix) >= self.maximum_feature:
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
