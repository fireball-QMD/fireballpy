#!/usr/bin/python3.8
import sys
import os
sys.path.append("/home/dani/fireballpy/build")
sys.path.append("/home/dani/fireballpy/geometry")
import fireballpy as fb

fb.loadfdata_from_path("/home/dani/Fdata_HC-new/")
fb.info_fdata()
fb.loadbas_from_file("/home/dani/fireballpy/test/input.bas")
fb.loadlvs_100()
fb.loadkpts_gamma()
fb.call_allocate_system()
fb.call_scf_loop()
print('  ========== CHARGES ====== ')
fb.info_charges()
fb.call_getenergy()
fb.info_energy()
fb.call_getforces()
fb.info_forces()



