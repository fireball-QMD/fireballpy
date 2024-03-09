here=$(pwd)
cd build
make
cd $here
cd test
../build/fireballpy.x
cd ..

