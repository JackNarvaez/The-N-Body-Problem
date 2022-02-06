/*---------------------------------------------
Evolution of a system of particles enclosed in
a box of unit sides, the only interaction is 
gravitation.
---------------------------------------------*/

#include <mpi.h>
#include <iostream>
#include <fstream> 
#include <sstream> 
#include <vector>
#include <random>
#include <cstdlib>
#include <cmath>
#include "NBodies.h"

int main(int argc, char **argv){

    int pId;  // Rank of process
    int nP;  // Number of processes
    int tag{0};
    int N; // Number of Particles
    int root{0};    // root process
    int steps = atoi(argv[1]);
    double dt = atof(argv[2]);
    double jump = atoi(argv[3]);
    MPI_Status status;
    std::string input = argv[4];
    read_NParticles(input, N);
    
    /*Initializes MPI*/
    MPI_Init(&argc, &argv);  
    MPI_Comm_size(MPI_COMM_WORLD, &nP);
    MPI_Comm_rank(MPI_COMM_WORLD, &pId);

    // Array size of each process
    std::vector<int> len(nP, 0);
    int end, begin;
    for(int ii=0; ii<nP; ii++){
        end = double(N)/nP*(ii+1);
        begin = double(N)/nP*ii;
        len[ii] = end-begin;
    }

    // Position, Velocity and Acceleration arrays
    std::vector<double> Pos(4*len[pId],0.0);    // [x, y, z, mass]
    std::vector<double> Vel(3*len[pId],0.0);    // [vx, vy, vz]
    std::vector<double> Acc(3*len[pId],0.0);    // [ax, ay, az]

    // Fill position and velocity vectors with the initial conditions
    Initial_state(input, Pos, Vel, len, N,  tag, pId, nP, root, status);
    // Calculate total acceleration felt by all particles
    Acceleration(Pos, Acc, len, N, tag, pId, nP, root, status);
    //Save the Position felt by all particles
    std::ofstream Data;
    Evolution(Data, Pos, Vel, Acc, len, N, tag, pId, nP, root, status, steps, dt, jump);

    /*Finalizes MPI*/
    MPI_Finalize();

    return 0;
}