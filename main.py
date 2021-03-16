from typing import List
import numpy as np

from breed import Breed
from features import FeatureHelper
from models import Genotype, Generation
from plot import PlotHelper



# initial data

parent1 = [
    [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
    [1, 0, 0, 1, 1, 0, 0, 1, 1, 0]
]

parent2 = [
    [0, 1, 1, 0, 0, 1, 1, 0, 0, 1],
    [0, 0, 0, 1, 0, 0, 0, 1, 0, 0]
]

parents = np.array([parent1, parent2])
# print(parents)
print(f"Initialized parent genotypes with: {parents.shape} shape")

possibilites = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
possibilites_for_gebv = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]

population = 3


def make_generation(parents: List):
    return Generation(
        index=0,
        genotypes=[Genotype(np.array(parent)) for parent in parents],
        population=2
    )


# PlotHelper.feature_individual_for_generation(
#     make_generation([parent1, parent2]),
#     FeatureHelper.gebv_feature
# )


breed = Breed(parents, possibilites, population, selection="gebv", max_age=1, max_population=36,
              possibilities_for_selection=possibilites_for_gebv)
last_generation_index = breed.evaluate(10)

PlotHelper.desirable_alleles_percent(breed.generations)
# print(last_generation_index)
# print(len(breed.generations[last_generation_index].genotypes))
# print(breed.generations[last_generation_index])

# PlotHelper.feature_individual_for_generation(
#     breed.generations[-1],
#     FeatureHelper.ohv_feature
# )
# PlotHelper.feature_individual_for_generation(
#     breed.generations[-1],
#     FeatureHelper.gebv_feature
# )

