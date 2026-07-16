if [ "$1" = "-n" ]; then
    make clean
fi
make -j
./build/ContinuumSystem params/extra_params.config