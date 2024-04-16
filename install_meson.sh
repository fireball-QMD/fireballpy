meson setup build --wipe
meson compile -C build
python -m numpy.f2py -m _fireball -c src/libf2py.f90 -Ibuild/src/libfireball.a.p -Lbuild/src -lfireball --link-lapack_opt --lower
mv _fireball*.so fireballpy
cd test
python test2.py
cd ..
