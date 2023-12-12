#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7

import sys
sys.path.append("../bin/f2py")
import f90

with open('fireball.in', 'r') as file:
    numeros = file.readline().split()

f90.variables.num1 = float(numeros[0])
f90.variables.num2 = float(numeros[1])
f90.suma()
print(f90.variables.resultado)


