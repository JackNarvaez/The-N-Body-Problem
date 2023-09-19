#!/bin/bash
#SBATCH -A ICT23_SMR3872
#SBATCH -p boost_usr_prod
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time 00:10:00
#SBATCH --mem=490000MB

#modulefiles to be loaded to have MPI on Leonardo

module purge
module load openmpi/4.1.4--gcc--11.3.0-cuda-11.8
module load cineca-hpyc

Nb=10000		# Number of bodies in the galaxy
Np=64           # Number of processes
i=0.1 		    # Inclination of the galaxy (rads/pi)
w=0.6 		    # Angle in the xy plane of the galaxy (rads/pi)
steps=10000 	# Evolution steps
jump=100		# Data storage interval
dt=0.001		# Time step
rad=5000		# Radius of Galaxy (AU)

mpirun -n ${Np} python3 galaxy.py ${Nb} ${i} ${w} ${rad}
mpic++ -std=c++17 Evolution.cpp NBodies.cpp -o Galaxy.x
mpirun -np ${Np} ./Galaxy.x ${steps} ${dt} ${jump} ${Nb}