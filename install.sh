here=$(pwd)
rm -fr build
mkdir build  
cd build
cmake ..
make
cd $here
cd test
../build/pyreball
cd ..

