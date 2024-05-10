rm -fr build
rm fireballpy/_fireball.*.so

if [ "$1" = "intel" ]; then
	if [ "$2" = "debug" ]; then
		CC=icx FC=ifx meson setup build -Dblas=mkl
	else
		CC=icx FC=ifx meson setup build -Dblas=mkl --buildtype custom
	fi
else
	if [ "$2" = "debug" ]; then
		meson setup build
	else
		meson setup build --buildtype custom
	fi
fi

meson compile -C build
meson install -C build

# Solo funcionara con 3.12
cd test
python3.12 test03.py
cd ..
