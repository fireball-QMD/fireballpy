rm -fr build

if [ "$1" = "intel" ]; then
	if [ "$2" = "debug" ]; then
		CC=icx FC=ifx meson setup build -Doptimization=3 -Dblas=mkl-dynamic-lp64-iomp -Dpython.install_env=auto
	else
		CC=icx FC=ifx meson setup build -Dblas=mkl-dynamic-lp64-iomp -Dpython.install_env=auto --buildtype custom
	fi
elif [ "$1" = "gnu" ]; then
	if [ "$2" = "debug" ]; then
		meson setup build -Doptimization=3 -Dpython.install_env=auto
	else
		meson setup build -Dpython.install_env=auto --buildtype custom
	fi
else
	echo "Must choose either 'intel' or 'gnu'"
	exit 1
fi

meson compile -C build
meson install -C build

# meson test -C build molecule periodic --verbose
