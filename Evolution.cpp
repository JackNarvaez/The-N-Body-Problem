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
    int tag{0};
    int root{0};    // root process
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
    int * len = (int *) malloc(nP * sizeof(int));
    // std::vector<int> len(nP, 0);
    int end, begin;
    for(int ii=0; ii<nP; ii++){
        end = double(N)/nP*(ii+1);
        begin = double(N)/nP*ii;
        len[ii] = end-begin;
    }
    // Position, Velocity and Acceleration arrays
    double * Pos = (double *) malloc(3*len[pId]*sizeof(double));    // [x, y, z]
    double * Vel = (double *) malloc(3*len[pId]*sizeof(double));    // [vx, vy, vz]
    double * Acc = (double *) malloc(3*len[pId]*sizeof(double));    // [ax, ay, az]
    double * Mass = (double *) malloc(len[pId]*sizeof(double));    // [m]
    
    for (int ii = 0; ii < 3*len[pId]; ii++) Acc[ii] = 0.0;

    // Fill position and velocity vectors with the initial conditions
    read_data(input, Pos, Vel, Mass);

    // Calculate total acceleration felt by all particles
    //Acceleration(Pos, Mass, Acc, len, N, tag, pId, nP, root, status);
    //Save the Position felt by all particles
    std::ofstream Data;
    Evolution(Data, Pos, Vel, Mass, Acc, len, N, tag, pId, nP, root, status, steps, dt, jump, Acceleration, PEFRL);
    /*Finalizes MPI*/
    MPI_Finalize();

    return 0;
}