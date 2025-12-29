from typing import TypeVar

T_co = TypeVar("T_co", covariant=True)
Extent = tuple[int | float, int | float, int | float, int | float]
