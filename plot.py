from typing import Callable, List

import numpy as np
import matplotlib.pyplot as plt

from models import Generation, Genotype, Feature
from features import FeatureHelper


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
    def gebv_plot(cls, generations: List[Generation], possibilities):
        xs = [generation.index for generation in generations]
        gebv_value = [[[FeatureHelper.gebv_feature(genotype, possibilities) for genotype in generation.genotypes], [generation.population]] for generation in generations]
        ys = [sum(i[0])/i[1][0] for i in gebv_value]
        ys_max = [max(i[0]) for i in gebv_value]
        ys_min = [min(i[0]) for i in gebv_value]
        #print([np.sum(np.array([max(genotype.matrix[0, i], genotype.matrix[1, i])*possibilities[i] for i in range(10)]) for genotype in generations[0].genotypes)])
        ys_border1 = []
        for generation in generations:
            a = []
            for i in range(10):
                b = []
                for genotype in generation.genotypes:
                    b.append(max([genotype.matrix[0, i], genotype.matrix[1, i]]))
                a.append(max(b))
            ys_border1.append(sum(a))

        '''
        histogram
        
        inf, sup = min(ys_min), max(ys_max)
        step = (sup - inf) / 20
        xs_histogram, ys_histogram = [], []
        ys_histogram = [(2*inf + step*(2*i+1)) / 2 for i in range(20)]
        for value in gebv_value:
        '''

        ys_border2 = []
        for generation in generations:
            a = []
            for i in range(10):
                b = []
                for genotype in generation.genotypes:
                    b.append(min([genotype.matrix[0, i], genotype.matrix[1, i]]))
                a.append(min(b))
            ys_border2.append(sum(a))

        plt.plot(xs, ys, c='r')
        plt.ylabel("GEBV")
        plt.xlabel("Generation number")
        plt.title("GS")
        plt.fill_between(xs, ys, ys_border1, color='white')
        plt.fill_between(xs, ys, ys_border2, color='white')
        plt.fill_between(xs, ys_min, ys_max, color='skyblue')
        plt.grid()
        plt.ylim([0, 12])
        ax = plt.axes()
        ax.set(facecolor='grey')
        plt.show()

    @classmethod
    def desirable_alleles_percent(cls, generations: List[Generation], selection):
        ys = [np.sum([genotype.matrix for genotype in generation.genotypes])/
              (np.prod(generation.genotypes[0].matrix.shape)*len(generation.genotypes)) for generation in generations]
        xs = [generation.index for generation in generations]

        plt.plot(xs, ys)
        plt.title(f'proportion of desirable alleles for {selection}')

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
