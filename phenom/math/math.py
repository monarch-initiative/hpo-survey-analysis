from typing import Set, Union, Sequence
from functools import reduce


# Union types
Num = Union[int, float]


def geometric_mean(values: Sequence[Num]) -> float:
    return reduce(lambda x,y: x*y, values)**(1/(len(values)))


def mean(values: Sequence[Num]) -> float:
    return reduce(lambda x,y: x+y, values)/(len(values))


def sum(values: Sequence[Num]) -> Num:
    return reduce(lambda x,y: x+y, values)