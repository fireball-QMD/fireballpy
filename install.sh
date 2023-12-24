here=$(pwd)
rm -fr build
rm -fr bin
mkdir build bin 
cd build
cmake ..
make
cd $here
cd test
./pyreball
cd ..

