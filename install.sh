here=$(pwd)
rm -fr build
mkdir build  
cd build
cmake ..
make
f2py3 -m fireballpy -c ../src/libf2py.f90 --fcompiler='ifort' -I. libfireballpy.a --link-lapack_opt
cd $here

cd test
./test.py
cd ..

