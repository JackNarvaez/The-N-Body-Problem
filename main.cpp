/*---------------------------------------------
Evolution of a system of particles enclosed in
a box of unit side, the unique interaction is 
gravitation.
---------------------------------------------*/

#include "mpi.h"
#include <iostream>
#include <fstream> 
#include <sstream> 
#include <vector>
#include <random>
#include <cstdlib>
#include <cmath>
#include "NBodies.h"

int main(int argc, char **argv){

    int pId;  // rank of process
    int nP;  //number of processes
    int tag{0};
    int N{200};
    int root{0};
    MPI_Status status;
    std::string input = "Data.txt";

    /*Initializes MPI*/
    MPI_Init(&argc, &argv);  
    MPI_Comm_size(MPI_COMM_WORLD, &nP);
    MPI_Comm_rank(MPI_COMM_WORLD, &pId);

    /*Array size of each process*/
    std::vector<int> len(nP, 0);
    int end, begin;
    for(int ii=0; ii<nP; ii++){
        end = double(N)/nP*(ii+1);
        begin = double(N)/nP*ii;
        len[ii] = end-begin;
    }

    /*Position, Momemtum and Force arrays*/
    std::vector<double> Pos(4*len[pId],0.0);
    std::vector<double> Mom(3*len[pId],0.0);
    std::vector<double> Force(3*len[pId],0.0);

    /*Fills position and momentum vectors with the initial conditions*/
    Initial_state(input, Pos, Mom, len, N,  tag, pId, nP, root, status);
    /*Calculates total force felt by all particles*/
    Total_Force(Pos, Force, len, N, tag, pId, nP, root, status);
    /*Prints the force felt by all particles*/
    std::ofstream Data;
    Save_data(Data, Force, len, N, tag, pId, nP, root, status);
    /*Finalizes MPI*/
    MPI_Finalize();

    return 0;
}