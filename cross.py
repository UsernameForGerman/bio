from typing import Union
import numpy as np

from .models import Animal, Bull, Cow, Genotype, Phenotype


def get_animal_cross(first: Animal, second: Animal) -> Union[Animal, Cow, Bull]:
    """
    Define new animal(genotype) created by crossbreeding first & second animals.
    :return:
    """

    new_genotype = Genotype(first.genotype.left, second.genotype.right)
    return Bull(new_genotype)


def get_phenotype(animal: Animal) -> Phenotype:
    """
    Define phenotype of animal by its genotype.
    """
    phenotype_value = np.mean(animal.genotype.left + animal.genotype.right)
    return Phenotype(phenotype_value)


