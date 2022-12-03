from typing import List
import numpy as np
from random import randint

from breed import Breed
from features import FeatureHelper
from models import Genotype, Generation
from plot import PlotHelper
import os
import csv


selection1 = 'gebv'
selection2 = 'ie534'

parents = np.array([list([randint(0, 1) for k in range(10)] for j in range(2)) for parent in range(30)])

# path = os.getcwd()
# with open(path+'\\fixtures\\bulls.csv') as file_genotypes:
#     reader = csv.reader(file_genotypes)
#     parents = []
#     a1, a2 = [], []
#     name = None
#     for i in reader:
#         if i[0] == name:
#             a1.append(int(i[1]))
#             a2.append(int(i[2]))
#         else:
#             if name != None:
#                 parents.append([a1, a2])
#                 a1, a2 = [], []
#             name = i[0]
#             a1.append(int(i[1]))
#             a2.append(int(i[2]))
#     print(parents)
#     parents = np.array(parents)
# print(f"Initialized parent genotypes with: {parents.shape} shape")
#
# with open(path+'\\fixtures\\parameters.txt', 'r') as file_parameters:
#     parameters = file_parameters.readlines()
#     parameteres_dict = {}
#     for line in parameters:
#         a = line.rsplit('=')
#         a = [value.rsplit()[0] for value in a]
#         parameteres_dict[a[0]] = a[1]
#     print(parameteres_dict)
#
possibilites = [0.5 for i in range(parents.shape[2])]  # r∈[0, 0.5]**(L-1)
possibilites_for_selection = [0.5 for i in range(parents.shape[2])]  # beta  -vector of weights

# breed = Breed(parents, possibilites, **parameteres_dict)
# last_generation_index = breed.evaluate(10)
# print(last_generation_index)
#
# PlotHelper.desirable_alleles_percent(breed.generations, selection1)
breed = Breed(parents, possibilites, population_of_progeny=2, selection="gebv", max_age=1, max_population=10000000000,
              possibilities_for_selection=possibilites_for_selection, crosses=15, puberty_period=0, self_pollinated=False)
last_generation_index = breed.evaluate(10)
PlotHelper.desirable_alleles_percent(breed.generations, selection2)
# PlotHelper.gebv_plot(breed.generations, possibilites_for_selection)
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

