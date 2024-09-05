from __future__ import annotations
from collections import deque

import numpy as np
from numpy.typing import NDArray


def alloc_integer(*size: int) -> NDArray[np.int64]:
    return np.zeros(size, dtype=np.int64, order="F")


def alloc_real(*size: int) -> NDArray[np.float64]:
    return np.zeros(size, dtype=np.float64, order="F")


def file_as_deque(fpath: str) -> deque[str]:
    with open(fpath, "r") as fp:
        raw = fp.read().splitlines()
    return deque(raw)


def read_integer_array(data: deque[str], width: int,
                       lines: int = 1) -> NDArray[np.int64]:
    result = np.empty((lines, width), dtype=np.int64)
    for i in range(lines):
        line = read_line(data)
        for j in range(width):
            result[i, j] = to_integer(line[j])
    return result if lines > 1 else result.reshape(-1,)


def read_real_array(data: deque[str], width: int,
                    lines: int = 1) -> NDArray[np.float64]:
    result = np.empty((lines, width), dtype=np.float64)
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


def to_integer(string: str) -> np.int64:
    return np.int64(string)


def to_real(string: str) -> np.float64:
    return np.float64(string.upper().replace("D", "E"))
