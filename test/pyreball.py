#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7

import os
import sys
import numpy as np
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)
import fortran_suma


# Leer los numeros desde el archivo fireball.in
with open('fireball.in', 'r') as file:
    num1, num2 = map(float, file.readline().split())
resultado=0
print(num1,num2)
fortran_suma.suma(num1, num2,resultado)
print("La suma es:", resultado)
