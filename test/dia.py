#!/usr/bin/env python

import sys
import os
sys.path.append("/home/dani/fireballpy/build")
import fireballpy as fb

fb.loadfdata_from_path("/home/dani/Fdata_HC-new/")
fb.info_fdata()

fb.loadbas_from_file("/home/dani/fireballpy/test/dia2.bas")
fb.loadlvs_from_file("/home/dani/fireballpy/test/dia2.lvs")
fb.loadkpts_from_file("/home/dani/fireballpy/test/dia2.kpts")
fb.set_icluster(0)
fb.set_gamma(0)
fb.rescal_structure(3.569)
#fb.print_atoms_positions()
fb.call_allocate_system()
fb.call_scf_loop()
fb.call_getenergy()

fb.info_charges()
fb.info_energy()
fb.info_forces()



