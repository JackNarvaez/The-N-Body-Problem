Np=$(nproc)
for ii in $(seq 1 ${Np}); do
    mpirun -np ${ii} --oversubscribe ./scaling.x "$1";
done > scaling.txt

