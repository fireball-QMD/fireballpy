from __future__ import annotations
from collections import deque

import numpy as np
from numpy.typing import NDArray

from ._types import real, integer


def alloc_integer(*size: int) -> NDArray[integer]:
    return np.zeros(size, dtype=integer, order="F")


def alloc_real(*size: int) -> NDArray[real]:
    return np.zeros(size, dtype=real, order="F")


def file_as_deque(fpath: str) -> deque[str]:
    with open(fpath, "r") as fp:
        raw = fp.read().splitlines()
    return deque(raw)


def read_integer_array(data: deque[str], width: int,
                       lines: int = 1) -> NDArray[integer]:
    result = np.empty((lines, width), dtype=integer)
    for i in range(lines):
        line = read_line(data)
        for j in range(width):
            result[i, j] = to_integer(line[j])
    return result if lines > 1 else result.reshape(-1,)


def read_real_array(data: deque[str], width: int,
                    lines: int = 1) -> NDArray[real]:
    result = np.empty((lines, width), dtype=real)
    for i in range(lines):
        line = read_line(data)
        for j in range(width):
            result[i, j] = to_real(line[j])
    return result if lines > 1 else result.reshape(-1,)


def read_line(data: deque[str]) -> list[str]:
    return data.popleft().split()


def skip_lines(data: deque[str], lines: int = 1) -> None:
    for _ in range(lines):
        _ = data.popleft()


def to_integer(string: str) -> integer:
    return integer(string)


def to_real(string: str) -> real:
    return real(string.upper().replace("D", "E"))
