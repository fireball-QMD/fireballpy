rm -fr build
meson setup build
cd build
ninja
cp libfireball.a.p/*.mod .
f2py3 -m _fireball -c ../src/libf2py.f90 -I. libfireball.a --link-lapack_opt
mv _fireball*.so ../fireballpy/
cd ../test
python3 test2.py
cd ..
