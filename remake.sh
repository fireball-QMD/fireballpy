here=$(pwd)
cd build
make
if command -v ifort 2>/dev/null; then
	f2py3 -m fireball -c ../src/libf2py.f90 --fcompiler='ifort' -I. libfireballpy.a --link-lapack_opt
else
	f2py3 -m fireball -c ../src/libf2py.f90 -I. libfireballpy.a --link-lapack_opt
fi
cd $here
cd test
#./dia.py
#./test.py
./test2.py # Creo que esta linea solo funciona para Dani asi que dont´t worry
cd ..
