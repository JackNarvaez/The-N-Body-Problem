#!/bin/bash
#SBATCH -A ICT23_SMR3872
#SBATCH -p boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH --time 00:10:00
#SBATCH --mem=490000MB

#modulefiles to be loaded to have MPI on Leonardo

module purge
module load openmpi/4.1.4--gcc--11.3.0-cuda-11.8
module load cineca-hpyc

Np=10		    # Number of processes
i=0.1 		  # Inclination of the galaxy (rads/pi)
w=0.6 		  # Angle in the xy plane of the galaxy (rads/pi)
steps=100 	# Evolution steps
dt=0.001		# Time step
rad=5000		# Radius of Galaxy (AU)

mpic++ -std=c++17 Evolution.cpp NBodies.cpp -o Weakscaling.x

rm weak.data

for p in 1 2 3 4; do
  Bd=$((10 ** ${p}))
  mpirun -np ${Np} python3 galaxy.py ${Bd} ${i} ${w} ${rad}
  mpirun -np ${Np} ./Weakscaling.x ${steps} ${dt} ${steps} ${Bd}
done >> weak.data
