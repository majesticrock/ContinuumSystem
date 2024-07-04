if [ "$1" = "-n" ]; then
    make clean
fi
make -j
n_mpi=1
n_omp=$(($(nproc)/(2*$n_mpi)))
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/ContinuumSystem params/params.config omega_D 0.2 5