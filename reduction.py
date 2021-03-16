from typing import List, Iterable
import math
import numpy as np
from models import Generation, Genotype, Feature
from features import FeatureHelper
from selection import Selection


class HerdReduction:
    def __init__(self, generation: Generation, max_population, selection: str = "gebv",
                 possibilities: List[float] = None):
        if selection == "pcv" and not possibilities:
            raise AttributeError("Не передан список вероятностей")
        self.generation = generation
        self.selection = selection
        self.freq = possibilities
        self.max_population = max_population

    def _ie534_reduction(self, crosses: int = 1) -> Generation:
        sorted_features = sorted(
            enumerate([FeatureHelper.gebv_feature(genotype) for genotype in self.generation.genotypes]),
            key=lambda value: value[1],
            reverse=True
        )  # sorts by gebv: integer
        parents = np.array(self.generation.genotypes)
        best_part_parents = parents[[value[0] for value in sorted_features[:math.ceil(len(sorted_features) * 0.05)]]]
        # choose fix number (len(self.generation.genotypes) * 0.05) of parents with the highest gebv
        len_best_part_parents = len(best_part_parents)
        parent_matrix = np.zeros((len_best_part_parents, len_best_part_parents))
        for i in range(len_best_part_parents):
            for j in range(len_best_part_parents):
                parent_matrix[i, j] = np.sum(
                    np.maximum(*np.maximum(best_part_parents[i].matrix, best_part_parents[j].matrix))
                )  # the number of locus where at least the one in pair have a desirable allele
                # * -  unpacking in arguments
        sorted_indexes = set()
        parent_matrix1 = parent_matrix.copy()
        while len(sorted_indexes) < self.max_population:
            index_parent1, index_parent2 = np.unravel_index(parent_matrix1.argmax(), parent_matrix1.shape)
            sorted_indexes.add(index_parent2)
            sorted_indexes.add(index_parent1)
            parent_matrix1[index_parent1][index_parent2] = -1
        if len(sorted_indexes) > self.max_population:
            sorted_indexes = np.array(sorted_indexes)
            v534 = np.array(sum([parent_matrix[i][j] for j in sorted_indexes if i != j]) for i in sorted_indexes)
            np.delete(sorted_indexes, np.argmin(v534))
        return Generation(
            index=self.generation.index,
            genotypes=[self.generation.genotypes[i] for i in sorted_indexes],
            population=self.max_population
        )

    def _gebv_reduction(self) -> Generation:
        sorted_features = sorted(
            enumerate([FeatureHelper.gebv_feature(genotype, self.freq) for genotype in self.generation.genotypes]),
            key=lambda value: value[1],
            reverse=True
        )
        indexes = [sorted_features[i][0] for i in range(self.max_population)]
        return Generation(
            index=self.generation.index,
            genotypes=[self.generation.genotypes[i] for i in indexes],
            population=self.max_population
        )

    def _ohv_reduction(self) -> Generation:
        sorted_features = sorted(
            enumerate([FeatureHelper.ohv_feature(genotype) for genotype in self.generation.genotypes]),
            key=lambda value: value[1],
        )
        indexes = [sorted_features[i][0] for i in range(self.max_population)]
        return Generation(
            index=self.generation.index,
            genotypes=[self.generation.genotypes[i] for i in indexes],
            population=self.max_population
        )

    def _pcv_selection(self) -> Generation:
        pcv_matrix = np.zeros((self.generation.population, self.generation.population))
        for i in range(self.generation.population):
            for j in range(i, self.generation.population):
                water_matrix = self.__water_matrix(self.generation.genotypes[i], self.generation.genotypes[j])
                pcv_matrix[i, j] = np.sum(water_matrix[:, :])
        sorted_indexes = set()
        pcv_matrix1 = pcv_matrix.copy()
        while len(sorted_indexes) < self.max_population:
            index_parent1, index_parent2 = np.unravel_index(pcv_matrix1.argmax(), pcv_matrix1.shape)
            sorted_indexes.add(index_parent2)
            sorted_indexes.add(index_parent1)
            pcv_matrix1[index_parent1][index_parent2] = -1
        if len(sorted_indexes) > self.max_population:
            sorted_indexes = np.array(sorted_indexes)
            pcv_sum = np.array(sum([pcv_matrix[i][j] for j in sorted_indexes if i != j]) for i in sorted_indexes)
            np.delete(sorted_indexes, np.argmin(pcv_sum))
        return Generation(
            index=self.generation.index,
            genotypes=[self.generation.genotypes[i] for i in sorted_indexes],
            population=self.max_population
        )

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

