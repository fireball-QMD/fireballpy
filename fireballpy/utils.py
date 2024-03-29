from collections import deque

from numpy.typing import NDArray
import numpy as np


def _alloc(dims: tuple[int, ...], dtype: str) -> NDArray:
    match dtype:
        case "float":
            return np.zeros(dims, dtype=float, order="F")
        case "int":
            return np.zeros(dims, dtype=int, order="F")
        case _:
            raise ValueError(
                f"dtype must be either 'int' or 'float'"
                f"got {dtype}"
            )


def _read_array(data: deque[str], dtype: str, lines: int) -> NDArray:
    val = data.popleft()
    match dtype:
        case "float":
            return np.array(
                [val.upper().replace("D", "E").split() for line in lines],
                dtype=float, order="F"
            )
        case "int":
            return np.array(
                [data.popleft().split() for line in lines],
                dtype=int, order="F"
            )
        case _:
            raise ValueError(
                f"dtype must be either 'int' or 'float'"
                f"got {dtype}"
            )


def _read_entry(data: deque[str], dtype: str, field: int) -> int | float | str:
    val = data.popleft()
    match dtype:
        case "float":
            return float(val.upper().replace("D", "E").split()[field])
        case "int":
            return int(val.split()[field])
        case "str":
            return val.split()[field]
        case _:
            raise ValueError(
                f"dtype must be either 'int', 'float' or 'str',"
                f"got {dtype}"
            )


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
    if len(dtypes) != len(line):
        raise ValueError(f"Wrong number of arguments, "
                         f"expected {len(line) + 1} got {len(dtypes) + 1}")
    for i, (l, t) in enumerate(zip(line, dtypes)):
        match t:
            case "float":
                line[i] = float(l.upper().replace("D", "E"))
            case "int":
                line[i] = int(l)
            case "str":
                continue
            case _:
                raise ValueError(
                    f"dtypes must be either 'int', 'float' or 'str'"
                    f", got {dtypes}")
    return line


def read_float_array(data: deque[str], lines: int = 1) -> NDArray[float]:
    return _read_array(data, "float", lines)


def read_int_array(data: deque[str], lines: int = 1) -> NDArray[int]:
    return _read_array(data, "int", lines)


def read_float_entry(data: deque[str], field: int = 0) -> float:
    return _read_entry(data, "float", field)


def read_int_entry(data: deque[str], field: int = 0) -> int:
    return _read_entry(data, "int", field)


def read_str_entry(data: deque[str], field: int = 0) -> str:
    return _read_entry(data, "str", field)


def skip_lines(data: deque[str], lines: int = 1):
    for _ in range(lines):
        _ = data.popleft()
