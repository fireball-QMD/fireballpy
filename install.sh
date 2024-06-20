rm -fr build
rm fireballpy/_fireball.*.so

if [ "$1" = "intel" ]; then
	if [ "$2" = "debug" ]; then
		CC=icx FC=ifx meson setup build -Dblas=mkl-dynamic-lp64-iomp
	else
		CC=icx FC=ifx meson setup build -Dblas=mkl-dynamic-lp64-iomp --buildtype custom
	fi
elif [ "$1" = "gnu" ]; then
	if [ "$2" = "debug" ]; then
		meson setup build
	else
		meson setup build --buildtype custom
	fi
else
	echo "Must choose either 'intel' or 'gnu'"
	exit 1
fi

meson compile -C build
meson install -C build

# meson test -C build molecule periodic --verbose
