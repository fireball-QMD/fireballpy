rm -fr build
meson setup build --buildtype debugoptimized
meson compile -C build
meson install -C build
cd test
python3.8 test02.py
cd ..
