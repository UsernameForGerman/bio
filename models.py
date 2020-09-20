# from __future__ import annotations
from typing import List
import pandas as pd


class Genotype:
    def __init__(self, left: List[int], right: List[int]):
        self.left = left
        self.right = right


class Phenotype:
    def __init__(self, value: int):
        self.value = value


class Animal:
    def __init__(self, name: str, genotype: Genotype, phenotype: Phenotype=None):
        self.genotype = genotype
        self.phenotype = phenotype
        self.name = name

    def __str__(self):
        return f'<{self.__class__.__name__}(name={self.name})>'

    def __repr__(self):
        return f'<{self.__class__.__name__}(name={self.name})>'

    @classmethod
    def read_csv(cls, filepath: str):
        genotype_values_df = pd.read_csv(filepath)
        animals = list()
        for group in genotype_values_df.groupby(genotype_values_df.bull):
            animals.append(
                cls(
                    name=group[0],
                    genotype=Genotype(
                        left=group[1].left,
                        right=group[1].right,
                    )
                )
            )
        return animals


class Bull(Animal):
    def __init__(self, name: str, genotype: Genotype, phenotype: Phenotype=None):
        super().__init__(name, genotype, phenotype)


class Cow(Animal):
    def __init__(self, name: str, genotype: Genotype, phenotype: Phenotype = None):
        super().__init__(name, genotype, phenotype)

