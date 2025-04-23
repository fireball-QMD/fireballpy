import numpy as np
from numpy.typing import NDArray

def generate_wavefunctions(ioption_in: np.int64, nexcite_in: np.int64, nznuc_in: np.int64, nzval_in: np.int64,
                           nzval_ion_in: np.int64, nzval_pp_in: np.int64, ioptim_in: np.int64,
                           atomname_in: str, ppfile_in: str, ppionfile_in: str, outpath_in: str,
                           sav_in: NDArray[np.bool_], lam_in: NDArray[np.int64], a0_in: NDArray[np.float64], rcutoff_in: NDArray[np.float64],
                           xocc_in: NDArray[np.float64], xocc0_in: NDArray[np.float64], xocc_ion_in: NDArray[np.float64],
                           cmix_in: NDArray[np.float64], r0_in: NDArray[np.float64], v0_in: NDArray[np.float64],
                           filename_wf_in: list[str], filename_ewf_in: list[str]) -> None: ...

def generate_vnn(nexcite_in: np.int64, nznuc_in: np.int64, ppfile_in: str, ppionfile_in: str, outpath_in: str,
                 lam_in: NDArray[np.int64], rcutoff_in: NDArray[np.float64], filename_wf_in: list[str], filename_ewf_in: list[str],
                 filename_na_in: list[str], filename_ena_in: list[str]) -> None: ...
