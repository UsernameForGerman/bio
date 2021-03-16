from typing import List, Callable
import numpy as np

from models import Generation, Genotype, Feature
from utils import filter_generation_by_index


class FeatureHelper:
    @classmethod
    def get_feature(cls, approach="gebv"):
        if approach == "gebv":
            return getattr(cls, "gebv_feature")

    @classmethod
    def gebv_feature(cls, genotype: Genotype, possibilities=None):
        return np.sum(genotype.matrix.dot(possibilities))

    @classmethod
    def ohv_feature(cls, genotype: Genotype):
        return np.sum(2*np.maximum(genotype.matrix[0], genotype.matrix[1]))

    @classmethod
    def pcv_feature(cls, genotype1: Genotype, genotype2: Genotype) -> int:
        return 0

    @classmethod
    def calculate_features_for_generations(cls, generations: List[Generation], approach: str = None,
                                           feature_function: Callable[[Genotype], Feature] = None) -> List[Feature]:
        _feature_function = cls.__get_valid_feature_function(approach, feature_function)
        generation_features = []
        for generation in generations:
            np.append(generation_features, cls._calculate_feature_for_generation(generation, _feature_function), axis=0)
        return generation_features

    @classmethod
    def calculate_feature_for_generation(
            cls, generations: List[Generation], generation_number: int = 0,
            approach: str = None, feature_function: Callable[[Genotype], Feature] = None) -> List[Feature]:
        _feature_function = cls.__get_valid_feature_function(approach, feature_function)
        generation = filter_generation_by_index(generations, generation_number)

        return cls._calculate_feature_for_generation(generation, _feature_function)

    @classmethod
    def _calculate_feature_for_generation(cls, generation: Generation, feature_function: Callable[[Genotype], Feature]
                                          ) -> List[Feature]:
        generation_features = []
        for genotype in generation.genotypes:
            generation_features.append(feature_function(genotype))
        return generation_features

    @classmethod
    def __get_valid_feature_function(cls, approach=None, feature_function=None):
        """
        Return selected feature function
        """
        if not approach and not feature_function:
            raise AttributeError("Выберите метод вычисления или укажите явно функцию вычисления")

        _feature_function = feature_function
        if not _feature_function:
            if not getattr(cls, f"{approach}_feature", None):
                raise NotImplementedError(f"Метод вычисления {approach} не реализован!")
            _feature_function = getattr(cls, f"{approach}_feature")

        return _feature_function
