if [ "$1" = "-n" ]; then
    make clean
fi
make -j
ITER_TYPE=g
ITER_MAX=0.2
ITER_NUM=10
./build/ContinuumSystem params/params.config $ITER_TYPE $ITER_MAX $ITER_NUM