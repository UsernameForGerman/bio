from typing import List, Iterable
import numpy as np

from models import Generation, Genotype
from features import FeatureHelper


class Selection:
    def __init__(self, generation: Generation, selection: str = "gebv", possibilites: List[float] = None):
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
        print(self.generation.population)
        print(self.generation.genotypes.__len__())
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

