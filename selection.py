from typing import List, Iterable
import numpy as np
import math
import itertools
from sys import maxsize


from models import Generation, Genotype
from features import FeatureHelper


class Selection:
    def __init__(self, generation: Generation, max_generations: int, selection: str = "gebv", possibilities: List[float] = None,
                 crosses: int = 1):
        if selection == "pcv" and not possibilities:
            raise AttributeError("Не передан список вероятностей")
        self.generation = generation
        self.selection = selection
        self.freq = possibilities
        self.crosses = crosses
        self.max_generations = max_generations

    def select(self) -> Iterable[Genotype]:
        return getattr(self, f"_{self.selection}_selection")(self.generation)
    
    def _ie534_selection(self, crosses: int = 1) -> List[Iterable[Genotype]]:
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
            for j in range(i): #excluding crossing with myself
                parent_matrix[i, j] = np.sum(
                    np.maximum(*np.maximum(best_part_parents[i].matrix, best_part_parents[j].matrix))
                )  # the number of locus where at least the one in pair have a desirable allele
                # * -  unpacking in arguments
        t = self.generation.index + 1
        K = 6 * (self.max_generations - t) # according to the model
        indexes = np.zeros(K)
        for i in range(K):
            index_parent1, index_parent2 = np.unravel_index(parent_matrix.argmax(), parent_matrix.shape)
            # i, j of the max
            indexes[i] = (index_parent1, index_parent2)
            parent_matrix[index_parent1, index_parent2] = -1
        return [(self.generation.genotypes[parent1], self.generation.genotypes[parent2]) for parent1, parent2 in indexes]

    def _gebv_selection(self) -> List[Iterable[Genotype]]:
        sorted_features = sorted(
            enumerate([FeatureHelper.gebv_feature(genotype, self.freq) for genotype in self.generation.genotypes]),
            key=lambda value: value[1],
            reverse=True
        )
        indexes = list()
        if len(sorted_features) // 2 < self.crosses:
            sorted_features = self.crosses_expansion(sorted_features)
        for i in range(self.crosses):
            indexes.append((sorted_features[2*i][0], sorted_features[2*i+1][0]))
        return [(self.generation.genotypes[parent1], self.generation.genotypes[parent2]) for parent1, parent2 in indexes]

    def _ohv_selection(self) -> List[Iterable[Genotype]]:
        sorted_features = sorted(
            enumerate([FeatureHelper.ohv_feature(genotype) for genotype in self.generation.genotypes]),
            key=lambda value: value[1],
        )
        indexes = list()
        if len(sorted_features) // 2 < self.crosses:
            sorted_features = self.crosses_expansion(sorted_features)
        for i in range(self.crosses):
            indexes.append((sorted_features[2*i][0], sorted_features[2*i+1][0]))
        return [(self.generation.genotypes[parent1], self.generation.genotypes[parent2]) for parent1, parent2 in indexes]

    def _pcv_selection(self): # -> List[Iterable[Genotype]]:

        '''
        pcv_matrix = np.zeros((self.generation.population, self.generation.population))
        for i in range(self.generation.population):
            for j in range(i, self.generation.population):
                water_matrix = self.__water_matrix(self.generation.genotypes[i], self.generation.genotypes[j])
                pcv_matrix[i, j] = np.sum(water_matrix[:, :])

        indexes = np.zeros(self.crosses)
        for i in range(self.crosses):
            index_parent1, index_parent2 = np.unravel_index(pcv_matrix.argmax(), pcv_matrix.shape)
            # i, j of the max
            indexes[i] = (index_parent1, index_parent2)
            pcv_matrix[index_parent1][index_parent2] = -1
        return [(self.generation.genotypes[parent1], self.generation.genotypes[parent2]) for parent1, parent2 in
                indexes]
        '''

        raise NotImplementedError('PCV not implemented')



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

    def crosses_expansion(self, sorted_features) -> tuple:
        difference = self.crosses - len(sorted_features)
        n = 2
        while math.factorial(n) / (
                2 * math.factorial(n - 2)) < difference:  # we select the required number of additional combinations
            n += 1
        addition = itertools.combinations(sorted_features[:n], 2)
        for i in addition:
            sorted_features.append(i[0])
            sorted_features.append(i[1])
        return sorted_features