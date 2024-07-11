if [ "$1" = "-n" ]; then
    make clean
fi
make -j
n_mpi=4
n_omp=$(($(nproc)/(2*$n_mpi)))
ITER_TYPE=phonon_coupling
ITER_MAX=39
ITER_NUM=8
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/ContinuumSystem params/params.config $ITER_TYPE $ITER_MAX $ITER_NUM