from typing import List, Iterable
import numpy as np
import math

from models import Generation, Genotype
from features import FeatureHelper


class Selection:
    def __init__(self, generation: Generation, selection: str = "gebv", possibilites: List[float] = None, max_number: int = sys.maxsize):
        if selection == "pcv" and not possibilites:
            raise AttributeError("Не передан список вероятностей")
        self.generation = generation
        self.selection = selection
        self.freq = possibilites
        self.max_number = max_number

    def select(self) -> Iterable[Genotype]:
        return getattr(self, f"_{self.selection}_selection")(self.generation)
    
    def _ie534_selection(self, crosses: int = 1) -> Iterable[Genotype]:
        sorted_features = sorted(
            enumerate([FeatureHelper.gebv_feature(genotype) for genotype in self.generation.genotypes]),
            key=lambda value: value[1],
            reverse=True
        )
        parents = np.array(self.generation.genotypes)
        best_part_parents = parents[[value[0] for value in sorted_features[:math.ceil(len(sorted_features) * 0.05)]]]
        len_best_part_parents = len(best_part_parents)
        parent_matrix = np.zeros((len_best_part_parents, len_best_part_parents))
        for i in range(len_best_part_parents):
            for j in range(len_best_part_parents):
                parent_matrix[i, j] = np.sum(
                    np.maximum(*np.maximum(best_part_parents[i].matrix, best_part_parents[j].matrix))
                )

        index_parent1, index_parent2 = np.unravel_index(parent_matrix.argmax(), parent_matrix.shape)
        return self.generation.genotypes[index_parent1], self.generation.genotypes[index_parent2]

    def _gebv_selection(self) -> Iterable[Genotype]:
        sorted_features = sorted(
            enumerate([FeatureHelper.gebv_feature(genotype, self.freq) for genotype in self.generation.genotypes]),
            key=lambda value: value[1],
            reverse=True
        )
        indexes = list()
        for i in range(0, self.max_number - 1, 2):
            indexes.append((sorted_features[i][0], sorted_features[i+1][0]))
        return [(self.generation.genotypes[i[0]], self.generation.genotypes[i[1]]) for i in indexes]

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
        n = self.generation.genotypes[0].matrix.shape[1]
        matrix = np.zeros((n, 4))
        parents = np.array([
            genotype1[0, :],
            genotype1[1, :],
            genotype2[0, :],
            genotype2[1, :]
        ])
        matrix[0, :] = parents[:, 0] / 4
        for i in range(1, matrix.shape[0]):
            t = self.__transition_matrix(i - 1)
            for j in range(matrix.shape[1]):
                matrix[i, j] = parents[j, i] * np.sum([t[k, j]*matrix[i - 1, k] for k in range(4)])

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

