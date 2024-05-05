rm -fr build
rm fireballpy/_fireball.cpython-38-x86_64-linux-gnu.so
meson setup build --buildtype debugoptimized
meson compile -C build
meson install -C build
cd test
python3.8 test03.py
cd ..
