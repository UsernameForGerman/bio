from typing import List
import numpy as np
from random import choices
import sys

from models import Genotype, Generation
from selection import Selection
from reduction import HerdReduction


class Breed:
    def __init__(self, parent_genotypes: np.ndarray, possibilites: List[float], population_of_progeny: int,
                 maximum_feature=None, selection: str = "gebv", max_age: int = 1, max_population: int = sys.maxsize,
                 possibilities_for_selection: List[float] = None, crosses: int = 1, puberty_period=0):
        if len(parent_genotypes.shape) != 3 or parent_genotypes.shape[1] != 2 or any(_ <= 0 for _ in parent_genotypes.shape):
            raise AttributeError("Массив генотипов особей задан неверно! Размерность должна быть (2 x 2 x N)")
        if max_age <= 0:
            raise AttributeError("Максимальный возраст сущности должен быть >= 1")
        if not Selection.selection_implemented(selection):
            raise NotImplementedError(f"Селекция {selection} не реализована!")
        if len(parent_genotypes) // 2 < crosses:
            raise AttributeError("Число скрещиваний не должно превышать количество родительских особей")
        if max_age <= puberty_period:
            raise AttributeError("Внимание! Ваша особь умирает прежде, чем начинает размножаться!")

        self.generations: List[Generation] = [Generation(
            index=0,
            genotypes=list(Genotype(genotype) for genotype in parent_genotypes),
            population=parent_genotypes.shape[0]
        )]
        self.possibilities = possibilites
        self.population_of_progeny = population_of_progeny
        self._current_generation = 0
        self.maximum_feature = (self.generations[0].genotypes[0].matrix.shape[1] * 2
                                if maximum_feature is None
                                else maximum_feature)
        self.selection = selection
        self.max_age = max_age
        self.max_population = max_population
        self.possibilities_for_selection = possibilities_for_selection
        self.crosses = crosses
        self.puberty_period = puberty_period

    def evaluate(self, max_generations: int = None):
        current_generation_number = 0
        while True:
            # if self._is_max_generation(current_generation_number):
            #     return current_generation_number

            parents_generation = self.filter_generation_for_puberty_period(self.get_generation(current_generation_number))
            child_generation = self.get_child_generation(parents_generation, max_generations)
            young_parents_genotypes = self.filter_generation_for_max_age(self.get_generation(current_generation_number))
            for parent in young_parents_genotypes:
                parent.age += 1
            child_generation.genotypes.extend(young_parents_genotypes)
            child_generation.population = len(child_generation.genotypes)
            if child_generation.population > self.max_population:
                reduction = HerdReduction(child_generation, self.max_population,
                                          possibilities=self.possibilities_for_selection)
                child_generation = getattr(reduction, f"_{self.selection}_selection")()
                child_generation.population = self.max_population
            self.generations.append(child_generation)

            if max_generations and current_generation_number == max_generations:
                return current_generation_number
            current_generation_number += 1

    def get_generation(self, generation_index: int) -> Generation:
        return list(filter(lambda generation: generation.index == generation_index, self.generations))[0]

    def get_child_generation(self, parent_generation: Generation, max_generations: int):
        if parent_generation.genotypes:
            selection = Selection(parent_generation, max_generations, possibilities=self.possibilities_for_selection, crosses=self.crosses)
            parents = getattr(selection, f"_{self.selection}_selection")()
            childrens = []
            for pair in parents:
                childrens += self.get_reproduce(pair[0], pair[1])
            #childrens = [self.get_reproduce(pair[0], pair[1]) for pair in parents]
        else:
            childrens = list()
        return Generation(index=parent_generation.index + 1, genotypes=childrens, population=len(childrens))

    def filter_generation_for_max_age(self, generation: Generation) -> List[Genotype]:
        return list(filter(lambda genotype: genotype.age < self.max_age, generation.genotypes))

    def get_reproduce(self, parent1: Genotype, parent2: Genotype) -> List[Genotype]:
        """
        Return children(which is amount=population_of_progeny)' genotype from parents' genotypes
        """
        children = []
        for _ in range(self.population_of_progeny):
            new_genotype = None
            for genotype in (parent1, parent2):
                if new_genotype is not None:
                    new_genotype = np.append(
                        new_genotype, [self.get_gamete(genotype)],
                        axis=0
                    )
                else:
                    new_genotype = np.array([self.get_gamete(genotype), ])
            children.append(Genotype(matrix=new_genotype))
        return children

    def get_genotype_indexies(self) -> List[int]:
        """
        Get indexies for choising haplotype in genotype
        """
        return_list = []
        for i in range(len(self.possibilities)):
            if i == 0:
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
            genotype[indexies[i], i] for i in range(genotype.matrix.shape[1])
        ])

    def _is_max_generation(self, generation_index: int):
        '''
        return True if at least one of the genotypes in generation is more than maximum_feature
        '''
        generation = self.get_generation(generation_index)
        for genotype in generation.genotypes:
            if self._is_max_feature_genotype(genotype):
                return True
        return False

    def _is_max_feature_genotype(self, genotype: Genotype):
        return True if np.sum(genotype.matrix) >= self.maximum_feature else False

    # def herd_reduction(self, generation: Generation) -> Generation:
    #
    #     pass
    #
    # def _is_max_number(self):
    #     number = 0
    #     for g in self.generations:
    #         number += g.population
    #     if number > self.max_number:
    #         self.herd_reduction()
    #
    def filter_generation_for_puberty_period(self, generation):
        new_generation = Generation(
            index = generation.index,
            genotypes = list(),
            population=0
        )
        for i in generation.genotypes:
            if i.age > self.puberty_period:
                new_generation.genotypes.append(i)
        return new_generation