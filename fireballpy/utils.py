from __future__ import annotations
from collections import deque

from numpy.typing import NDArray
import numpy as np


def _alloc(dims: tuple[int, ...], dtype: str) -> NDArray:
    if dtype == "float":
        return np.zeros(dims, dtype=float, order="F")
    elif dtype == "int":
        return np.zeros(dims, dtype=int, order="F")
    else:
        raise ValueError(f"dtype must be either 'int' or 'float' got {dtype}")


def _read_array(data: deque[str],
                dtype: str,
                lines: int,
                width: None | int) -> NDArray:
    val = data.popleft()
    if dtype == "float":
        arr = np.array([val.upper().replace("D", "E").split()
                        for _ in range(lines)],
                       dtype=float, order="F")
    elif dtype == "int":
        arr = np.array([val.split() for _ in range(lines)],
                       dtype=int, order="F")
    else:
        raise ValueError(f"dtype must be either 'int' or 'float' got {dtype}")
    width = arr.shape[1] if width is None else width
    return arr[:width]


def _read_entry(data: deque[str], dtype: str, field: int) -> int | float | str:
    val = data.popleft()
    if dtype == "float":
        return float(val.upper().replace("D", "E").split()[field])
    elif dtype == "int":
        return int(val.split()[field])
    elif dtype == "str":
        return val.split()[field]
    else:
        raise ValueError(f"dtype must be either 'int', 'float' "
                         f"or 'str' got {dtype}")


def alloc_float(*size: int) -> NDArray[float]:
    return _alloc(size, "float")


def alloc_int(*size: int) -> NDArray[int]:
    return _alloc(size, "int")


def read_file(fpath: str, header: int = 0) -> deque[str]:
    with open(fpath, "r") as fp:
        raw = fp.read().splitlines()[header:]
    return deque(raw)


def read_line(data: deque[str], *dtypes: str) -> list[int | float | str]:
    line = data.popleft().split()
    if not dtypes:
        dtypes = tuple("str" for _ in line)
    if len(dtypes) > len(line):
        raise ValueError("Wrong number of arguments, "
                         f"{len(dtypes) + 1} < {len(line) + 1}")
    for i, (l, t) in enumerate(zip(line, dtypes)):
        if t == "float":
            line[i] = float(l.upper().replace("D", "E"))
        elif t == "int":
            line[i] = int(l)
        elif t == "str":
            continue
        else:
            raise ValueError(f"dtypes must be either 'int', 'float' "
                             f"or 'str' got {t} at position {i}")
    return line


def read_float_array(data: deque[str],
                     lines: int = 1,
                     width: None | int = None) -> NDArray[float]:
    return _read_array(data, "float", lines, width)


def read_int_array(data: deque[str],
                   lines: int = 1,
                   width: None | int = None) -> NDArray[int]:
    return _read_array(data, "int", lines, width)


def read_float_entry(data: deque[str], field: int = 0) -> float:
    return _read_entry(data, "float", field)


def read_int_entry(data: deque[str], field: int = 0) -> int:
    return _read_entry(data, "int", field)


def read_str_entry(data: deque[str], field: int = 0) -> str:
    return _read_entry(data, "str", field)


def skip_lines(data: deque[str], lines: int = 1):
    for _ in range(lines):
        _ = data.popleft()
