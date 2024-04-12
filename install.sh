here=$(pwd)
rm -fr build
mkdir build
cd build
cmake ..
make
if command -v ifort 2>/dev/null; then
	f2py3 -m _fireball -c ../src/libf2py.f90 --fcompiler='ifort' -I. libfireballpy.a --link-lapack_opt
else
	f2py3 -m _fireball -c ../src/libf2py.f90 -I. libfireballpy.a --link-lapack_opt
fi
mv _fireball*.so ../fireballpy/
cd $here

cd test
python3 test2.py
cd ..
