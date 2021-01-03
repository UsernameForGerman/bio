from typing import Callable, List

import numpy as np
import matplotlib.pyplot as plt

from models import Generation, Genotype, Feature


class PlotHelper:

    @classmethod
    def _feature_individual(cls, genotype: Genotype, feature_function: Callable[[Genotype], int], index: int=0):
        return Feature(
            index=index,
            number=feature_function(genotype)
        )

    @classmethod
    def feature_individual_for_generation(cls, generation: Generation, feature_function: Callable[[Genotype], int]):
        """
        Plot graphic
        """
        features = []
        i = 0
        for genotype in generation.genotypes:
            features.append(cls._feature_individual(genotype, feature_function, i))
            i += 1

        plt.bar(
            [feature.index for feature in features],
            [feature.number for feature in features],
            width=0.3
        )
        plt.ylabel(feature_function.__name__)
        plt.xlabel("Характеристическое число")
        plt.title(f"Характеристическое число для поколения №{generation.index}")
        plt.show()

        plt.close()

    @classmethod
    def desirable_alleles_percent(cls, generations: List[Generation]):
        ys = [np.sum([genotype.matrix for genotype in generation.genotypes])/(np.prod(generation.genotypes[0].matrix.shape)*len(generation.genotypes)) for generation in generations]
        xs = [generation.index for generation in generations]

        plt.plot(xs, ys)
        plt.show()

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
