rm -fr variables.mod
echo "compilamos con ifort"
ifort -c ../src/variables.f90
ifort -c ../src/diagonalize_matrix.f90
ifort -c ../src/pyreball.f90 
ifort -o test.x variables.o diagonalize_matrix.o pyreball.o -mkl
rm -fr *.o
rm -fr *.mod
./test.x
echo "python"
python3 dia.py
rm -fr test.x


