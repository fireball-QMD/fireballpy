#!/usr/bin/env python

import sys
import os
sys.path.append("/home/dani/fireballpy/build")
import fireballpy

fireballpy.info()
fireballpy.loadfdata_from_file("/home/dani/Fdata_HC-new/")
fireballpy.loadbasformat_from_file("/home/dani/fireballpy/test/input.bas")
