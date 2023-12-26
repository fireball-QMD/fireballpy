import numpy as np

def read_complex_matrix(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        n = len(lines)
        matrix = np.zeros((n, n), dtype=np.complex128)
        for i, line in enumerate(lines):
            a=line.replace('(','').replace(')','').replace(',','').split()
            #print(a[0],a[1])
            for j in range(4):
              real=a[0+j*2]
              imag=a[1+j*2]
              matrix[i, j] = complex(float(real), float(imag))
    return matrix


file_name = 'fireball.in'
complex_matrix = read_complex_matrix(file_name)
eigenvalues, eigenvectors = np.linalg.eig(complex_matrix)
print(complex_matrix)
print("Autovalores:")
for eigenvalue in eigenvalues:
    print(f"{eigenvalue:.3f}")
