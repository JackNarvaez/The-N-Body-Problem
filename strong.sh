#!/bin/bash
#SBATCH -A ICT23_SMR3872
#SBATCH -p boost_usr_prod
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time 00:10:00
#SBATCH --mem=490000MB

# module files to be loaded to have MPI on Leonardo

module purge
module load openmpi/4.1.4--gcc--11.3.0-cuda-11.8
module load cineca-hpyc

Ng=10000		# Number of bodies in the galaxy
i=0.1 		  # Inclination of the galaxy (rads/pi)
w=0.6 		  # Angle in the xy plane of the galaxy (rads/pi)
steps=100 	# Evolution steps
dt=0.001		# Time step
rad=5000		# Radius of Galaxy (AU)

rm strong.data

mpic++ -std=c++17 Evolution.cpp NBodies.cpp -o Strongscaling.x

for Np in 1 2 4 8 16 32 64; do
  mpirun -np ${Np} python3 galaxy.py ${Ng} ${i} ${w} ${rad}
  mpirun -np ${Np} ./Strongscaling.x ${steps} ${dt} ${steps} ${Ng}
done >> strong.data