from typing import List

from models import Generation


def filter_generation_by_index(generations: List[Generation], generation_index: int = 0):
    generation_list = list(filter(lambda generation: generation.index == generation_index, generations))
    if not generation_list:
        raise AttributeError(f"Поколения №{generation_index} не существует в переданном массиве")
    return generation_list[0]











