#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7


import sys
import os
sys.path.append("/home/dani/fireballpy/build")
import fireballpy as fb

fb.loadfdata_from_path("/home/dani/Fdata_HC-new/")
fb.info_fdata()
fb.loadbas_from_file("/home/dani/fireballpy/test/input.bas")
fb.loadlvs_100()
fb.loadkpts_gamma()
fb.call_allocate_system()
fb.call_scf_loop()
fb.call_getenergy()
fb.info_energy()
fb.info_forces()

