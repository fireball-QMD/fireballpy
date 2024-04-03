here=$(pwd)
cd build
make
f2py3 -m fireballpy -c ../src/libf2py.f90 --fcompiler='ifort' -I. libfireballpy.a --link-lapack_opt
cd $here
cd test
#./dia.py
./test.py
cd ..

