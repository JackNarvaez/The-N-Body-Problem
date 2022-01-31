for N  in $(seq 1 200 8000); do
    g++ random.cpp -o random.x
    ./random.x ${N}
    mpirun -np "$1" ./scaling.x 1;
done > strong"$1".txt

