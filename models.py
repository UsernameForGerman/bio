from typing import List, Any
from nptyping import NDArray
from dataclasses import dataclass, fields


@dataclass
class DefaultVal:
    val: Any


@dataclass
class DefaultValuesMixin:
    def __post_init__(self):
        for field in fields(self):
            if isinstance(field.default, DefaultVal):
                field_val = getattr(self, field.name)
                if isinstance(field_val, DefaultVal) or field_val is None:
                    setattr(self, field.name, field.default.val)


@dataclass
class Genotype:
    matrix: NDArray
    age: int = 0

    def __getitem__(self, item):
        return self.matrix.__getitem__(item)


@dataclass
class Generation(DefaultValuesMixin):
    index: int
    genotypes: List[Genotype]
    population: int = 1


@dataclass
class Feature:
    index: int
    number: int
