for N  in $(seq 1 10 1000); do
    g++ random.cpp -o random.x
    ./random.x ${N}
    mpirun -np 1 ./scaling.x 1;
done > strong.txt

