from numpy import str_, int32, float64
from numpy.typing import NDArray

def generate_wavefunctions(nznuc_in: int32, nexcite_in: int32, nzval_in: int32,
                           nzval_ion_in: int32, ioptim_in: int32, nzval_pp_in: float64,
                           ppfile_in: str, ppionfile_in: str, outpath_in: str, sav_in: list[bool],
                           lam_in: NDArray[int32], rcutoff_in: NDArray[float64], xocc_in: NDArray[float64],
                           xocc0_in: NDArray[float64], xocc_ion_in: NDArray[float64], cmix_in: NDArray[float64],
                           exc_names_in: NDArray[str_], exc_weights_in: NDArray[float64],
                           r0_in: NDArray[float64], v0_in: NDArray[float64]) -> None: ...

def generate_vnn(nexcite_in: int32, nznuc_in: int32, ppfile_in: str,
                 ppionfile_in: str, outpath_in: str,
                 lam_in: NDArray[int32], rcutoff_in: NDArray[float64]) -> None: ...
