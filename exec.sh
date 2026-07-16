if [ "$1" = "-n" ]; then
    make clean
fi
make -j
./build/default/ContinuumSystem params/params.config