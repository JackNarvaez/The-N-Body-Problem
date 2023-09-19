/*---------------------------------------------
Evolution of the Sagittarius A* cluster.
Interactions are only due to gravitation.
.--------------------------------------------*/

#include <mpi.h>
#include <iostream>
#include <fstream> 
#include <sstream> 
#include <cstdlib>
#include <cmath>
#include <string>
#include "NBodies.h"

int main(int argc, char **argv){
   
    int pId;                        // Process rank
    int nP;                         // Number of processes
    int tag{0};                     // Tag message
    int root{0};                    // Root process
    int steps = atoi(argv[1]);      // Evolution steps
    double dt = atof(argv[2]);      // Time step
    int jump  = atoi(argv[3]);      // Data storage interval
    int N     = 14;                 // Total bodies
    MPI_Status status;
    body bd;                        // Bodies

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
    bd.m = (double *) malloc(len[pId]*sizeof(double));      // [m]

    for (int ii = 0; ii < 3*len[pId]; ii++) bd.a[ii] = 0.0;

    // Fill position and velocity vectors with the initial conditions
    if (pId == root) {
        std::string input = "SagA.data";
        body bd_glb;
        bd_glb.r = (double *) malloc(3*N*sizeof(double));    // [x, y, z]
        bd_glb.v = (double *) malloc(3*N*sizeof(double));    // [vx, vy, vz]
        bd_glb.m = (double *) malloc(N*sizeof(double));      // [m]

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

    std::ofstream Data;
    Evolution(Data, bd.r, bd.v, bd.m, bd.a, len, N, tag, pId, nP, root, status, steps, dt, jump, Acceleration, PEFRL);
    
    /*Finalizes MPI*/
    MPI_Finalize();

    return 0;
}