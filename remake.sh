here=$(pwd)
cd build
make
cd $here
cd test
../build/pyreball
cd ..

