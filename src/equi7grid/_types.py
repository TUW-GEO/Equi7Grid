# Copyright (c) 2026, TU Wien
# Licensed under the MIT License. See LICENSE file.

from typing import TypeVar

T_co = TypeVar("T_co", covariant=True)
Extent = tuple[int | float, int | float, int | float, int | float]
