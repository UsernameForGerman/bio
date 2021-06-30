from typing import List
import numpy as np
from random import randint

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

selection = 'gebv'

parents = np.array([list([randint(0, 1) for k in range(10)] for j in range(2)) for parent in range(30)])
print(parents)
print(f"Initialized parent genotypes with: {parents.shape} shape")

possibilites = [0.5 for i in range(10)]  # r∈[0, 0.5]**(L-1)
possibilites_for_selection = [0.5 for i in range(10)]  # beta  -vector of weights

population = 20


# def make_generation(parents: List):
#     return Generation(
#         index=0,
#         genotypes=[Genotype(np.array(parent)) for parent in parents],
#         population=2
#     )


# PlotHelper.feature_individual_for_generation(
#     make_generation([parent1, parent2]),
#     FeatureHelper.gebv_feature
# )


breed = Breed(parents, possibilites, population, selection="gebv", max_age=1, max_population=1000000000,
              possibilities_for_selection=possibilites_for_selection, crosses=15, puberty_period=0)
last_generation_index = breed.evaluate(10)
print(last_generation_index)

#PlotHelper.desirable_alleles_percent(breed.generations, selection)
PlotHelper.gebv_plot(breed.generations, possibilites_for_selection)
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

