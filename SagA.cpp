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
#include <chrono>
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
    body bd;
    int N = 14; // Number of total particles
    /*Initializes MPI*/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nP);
    MPI_Comm_rank(MPI_COMM_WORLD, &pId);

    // Length of each proccess' vectors
    int * len = (int *) malloc(nP * sizeof(int));
    len[0] = 7;
    len[1] = 7;

    // Position, Velocity and Acceleration arrays
    bd.r = (double *) malloc(3*len[pId]*sizeof(double));    // [x, y, z]
    bd.v = (double *) malloc(3*len[pId]*sizeof(double));    // [vx, vy, vz]
    bd.a = (double *) malloc(3*len[pId]*sizeof(double));    // [vx, vy, vz]
    bd.m = (double *) malloc(len[pId]*sizeof(double));    // [m]

    for (int ii = 0; ii < 3*len[pId]; ii++) bd.a[ii] = 0.0;

    if (pId == root) {
        std::string input = "SagA.data";
        body bd_glb;
        bd_glb.r = (double *) malloc(3*N*sizeof(double));    // [x, y, z]
        bd_glb.v = (double *) malloc(3*N*sizeof(double));    // [vx, vy, vz]
        bd_glb.m = (double *) malloc(N*sizeof(double));    // [m]

        read_data(input, bd_glb.r, bd_glb.v, bd_glb.m);

        MPI_Send(bd_glb.r, 3*len[pId], MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);
        MPI_Send(bd_glb.v, 3*len[pId], MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);
        MPI_Send(bd_glb.m, len[pId], MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);

        bd.r = &bd_glb.r[3*len[pId]];
        bd.v = &bd_glb.v[3*len[pId]];
        bd.m = &bd_glb.m[len[pId]];

    } else {
        MPI_Recv(bd.r, 3*len[pId], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(bd.v, 3*len[pId], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(bd.m, len[pId], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    }
    // Fill position and velocity vectors with the initial conditions

    // Calculate total acceleration felt by all particles
    //Acceleration(Pos, Mass, Acc, len, N, tag, pId, nP, root, status);
    //Save the Position felt by all particles
    std::ofstream Data;
    Evolution(Data, bd.r, bd.v, bd.m, bd.a, len, N, tag, pId, nP, root, status, steps, dt, jump, Acceleration, PEFRL);
    /*Finalizes MPI*/

    MPI_Finalize();

    return 0;
}