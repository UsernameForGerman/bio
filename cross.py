from typing import Union, List, Iterable, Dict, Callable, Any
import numpy as np
from nptyping import NDArray
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
    matrix: NDArray

    def __getitem__(self, item):
        return self.matrix.__getitem__(item)


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
            return getattr(cls, "gebv_feature")

    @classmethod
    def gebv_feature(cls, genotype: Genotype):
        return np.sum(genotype.matrix)

    @classmethod
    def ohv_feature(cls, genotype: Genotype):
        return np.sum(2*np.maximum(genotype.matrix[0], genotype.matrix[1]))

    @classmethod
    def pcv_feature(cls, genotype1: Genotype, genotype2: Genotype) -> int:
        return 0

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
        generation = filter_generation_by_index(generations, generation_number)

        return cls._calculate_feature_for_generation(generation, _feature_function)

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


class Selection:
    def __init__(self, generation: Generation, selection: str = "gebv", possibilites: List[int] = None):
        if selection == "pcv" and not possibilites:
            raise AttributeError("Не передан список вероятностей")
        self.generation = generation
        self.selection = selection
        self.freq = possibilites

    def select(self) -> Iterable[Genotype]:
        return getattr(self, f"_{self.selection}_selection")(self.generation)

    def _gebv_selection(self) -> Iterable[Genotype]:
        sorted_features = sorted(
            enumerate([FeatureHelper.gebv_feature(genotype) for genotype in self.generation.genotypes]),
            key=lambda value: value[1],
        )
        index_parent1, index_parent2 = sorted_features[-1][0], sorted_features[-2][0]
        return self.generation.genotypes[index_parent1], self.generation.genotypes[index_parent2]

    def _ohv_selection(self) -> Iterable[Genotype]:
        sorted_features = sorted(
            enumerate([FeatureHelper.ohv_feature(genotype) for genotype in self.generation.genotypes]),
            key=lambda value: value[1],
        )
        index_parent1, index_parent2 = sorted_features[-1][0], sorted_features[-2][0]
        return self.generation.genotypes[index_parent1], self.generation.genotypes[index_parent2]

    def _pcv_selection(self) -> Iterable[Genotype]:
        pcv_matrix = np.zeros((self.generation.population, self.generation.population))
        for i in range(self.generation.population):
            for j in range(i, self.generation.population):
                water_matrix = self.__water_matrix(self.generation.genotypes[i], self.generation.genotypes[j])
                pcv_matrix[i, j] = np.sum(water_matrix[:, :])

        index_parent1, index_parent2 = np.unravel_index(pcv_matrix.argmax(), pcv_matrix.shape)
        return self.generation.genotypes[index_parent1], self.generation.genotypes[index_parent2]

    def __water_matrix(self, genotype1: Genotype, genotype2: Genotype):
        n = self.generation.genotypes[0].matrix.shape[0]
        matrix = np.zeros((n, 4))
        parents = np.concatenate((
            genotype1[:, 0],
            genotype1[:, 1],
            genotype2[:, 0],
            genotype2[:, 1]
        ))

        matrix[0, :] = 0.25 * parents[0, :]
        for i in range(1, matrix.shape[0]):
            t = self.__transition_matrix(i - 1)
            for j in range(matrix.shape[1]):
                matrix[i, j] = parents[i, j] * np.sum([t[k, j, i - 1]*matrix[i - 1, k] for k in range(4)])

        return matrix

    def __transition_matrix(self, i: int):
        """
        i - index of possibilities
        """
        return np.array([
            [(1 - self.freq[i])**2, self.freq[i]*(1 - self.freq[i]), 0.5*self.freq[i], 0.5*self.freq[i]],
            [self.freq[i]*(1 - self.freq[i]), (1 - self.freq[i])**2, 0.5*self.freq[i], 0.5*self.freq[i]],
            [0.5*self.freq[i], 0.5*self.freq[i], (1 - self.freq[i])**2, self.freq[i]*(1 - self.freq[i])],
            [0.5*self.freq[i], 0.5*self.freq[i], (1 - self.freq[i])**2, self.freq[i]*(1 - self.freq[i])]
        ])

    @classmethod
    def selection_implemented(cls, selection: str):
        return True if f"_{selection}_selection" in dir(cls) else False


class Breed:
    def __init__(self, parent_genotypes: np.ndarray, possibilites: List[float], population_of_progeny: int,
                 maximum_feature=None, selection: str = "gebv"):
        if len(parent_genotypes.shape) != 3 or parent_genotypes.shape[0:] != (2, 2) or parent_genotypes.shape[0] <= 0:
            raise AttributeError("Массив генотипов особей задан неверно! Размерность должна быть (2 x 2 x N)")
        if not Selection.selection_implemented(selection):
            raise NotImplementedError(f"Селекция {selection} не реализована!")

        self.generations: List[Generation] = [Generation(
            index=0,
            genotypes=list(Genotype(genotype) for genotype in parent_genotypes),
            population=parent_genotypes.shape[2]
        )]
        self.possibilities = possibilites
        self.population_of_progeny = population_of_progeny
        self._current_generation = 0
        self.maximum_feature = (self.generations[0].genotypes[0].matrix.shape[0] * 2
                                if maximum_feature is None
                                else maximum_feature)
        self.selection = selection

    def evaluate(self, max_generations: int = None):
        current_generation_number = 0
        while True:
            if self._is_max_generation(current_generation_number):
                return current_generation_number

            parents = self.get_generation(current_generation_number)
            children = self.get_child_generation(parents)
            self.generations.append(children)

            if max_generations and current_generation_number == max_generations:
                return current_generation_number
            current_generation_number += 1

    def get_generation(self, generation_index: int) -> Generation:
        return list(filter(lambda generation: generation.index == generation_index, self.generations))[0]

    def get_child_generation(self, parent_generation: Generation):
        parent1, parent2 = getattr(FeatureHelper, f"_{self.selection}_selection")(parent_generation)
        children = self.get_reproduce(parent1, parent2)
        return Generation(index=parent_generation.index + 1, genotypes=children, population=len(children))

    def get_reproduce(self, parent1: Genotype, parent2: Genotype) -> List[Genotype]:
        """
        Return children(which is amount=population_of_progeny)' genotype from parents' genotypes
        """
        children = []
        for _ in range(self.population_of_progeny):
            new_genotype = None
            for genotype in (parent1, parent2):
                if new_genotype:
                    np.append(
                        new_genotype, self.get_gamete(genotype),
                        axis=0
                    )
                else:
                    new_genotype = self.get_gamete(genotype)
            children.append(Genotype(matrix=new_genotype))
        return children

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

    def get_gamete(self, genotype: Genotype) -> np.ndarray:
        """
        Return gamete from genotype input
        """
        indexies = self.get_genotype_indexies()
        return np.array([
            genotype[i, indexies[i]] for i in range(genotype.matrix.shape[0])
        ])

    def _is_max_generation(self, generation_index: int):
        generation = self.get_generation(generation_index)
        for genotype in generation.genotypes:
            if self._is_max_feature_genotype(genotype):
                return True
        return False

    def _is_max_feature_genotype(self, genotype: Genotype):
        return True if np.sum(genotype.matrix) >= self.maximum_feature else False


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
