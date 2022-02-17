/*---------------------------------------------
Evolution of a system of particles, the only 
interaction is gravitation.
.--------------------------------------------*/

#include <mpi.h>
#include <iostream>
#include <fstream> 
#include <sstream> 
#include <vector>
#include <random>
#include <cstdlib>
#include <cmath>
#include <string>
#include "NBodies.h"

int main(int argc, char **argv){

    int pId;  // Rank of process
    int nP;  // Number of processes
    int tag = 0;
    int root = 0;    // root process
    int steps = atoi(argv[1]);
    double dt = atof(argv[2]);
    int jump = atoi(argv[3]);
    MPI_Status status;
    int N = atoi(argv[4]); // Number of total particles
    /*Initializes MPI*/
    MPI_Init(&argc, &argv);  
    MPI_Comm_size(MPI_COMM_WORLD, &nP);
    MPI_Comm_rank(MPI_COMM_WORLD, &pId);
    /*Reads number of particles in each file*/
    std::string data = "data";
    std::string id = std::to_string(pId);
    std::string extension = ".txt";
    std::string input = data + id + extension;
    // Length of each proccess' vectors
    std::vector<int> len(nP, 0);
    int end, begin;
    for(int ii=0; ii<nP; ii++){
        end = double(N)/nP*(ii+1);
        begin = double(N)/nP*ii;
        len[ii] = end-begin;
    }
    // Position, Velocity and Acceleration arrays
    std::vector<double> Pos;    // [x, y, z, mass]
    std::vector<double> Vel;    // [vx, vy, vz]
    std::vector<double> Acc(3*len[pId],0.0);   // [ax, ay, az]
    // Fill position and velocity vectors with the initial conditions
    read_data(input, Pos, Vel);
    // Calculate total acceleration felt by all particles
    Acceleration(Pos, Acc, len, N, tag, pId, nP, root, status);
    //Save the Position felt by all particles
    std::ofstream Data;
    Evolution(Data, Pos, Vel, Acc, len, N, tag, pId, nP, root, status, steps, dt, jump, Acceleration);
    /*Finalizes MPI*/
    MPI_Finalize();

    return 0;
}
