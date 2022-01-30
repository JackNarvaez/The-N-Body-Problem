for nP in $(seq 1 16); do
    mpirun -np ${nP} --oversubscribe ./scaling.x "$1";
done > scaling.txt

