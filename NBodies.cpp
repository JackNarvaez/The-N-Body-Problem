#include "mpi.h"
#include <iostream>
#include <vector>
#include <random>
#include <cstdlib>
#include <cmath> 

void Initial_configuration(std::vector<double>&Pos, std::vector<double>&Mom, std::vector<int>&len, int N, int tag, int id, int np, MPI_Status status);
void Total_Force(std::vector<double>& Pos, std::vector<double>& Force, std::vector<int>& len, int N, int tag, int id, int np, MPI_Status status);
void Gravitational_force(std::vector<double>& Force, const std::vector<double>& vec0, const std::vector<double>& vec1, int len0, int len1);
void print_all(std::vector<double>& Vec, const std::vector<int>& len, int N,  int tag, int id, int np, MPI_Status status);
void print_vec(std::vector<double>& Vec, int len);

int main(int argc, char **argv){
  /*Variables*/
  int id, np, tag=0, N=20;
  MPI_Status status;
  /*Initializes MPI*/
  MPI_Init(&argc, &argv);  
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  /*Size of each process array*/
  std::vector<int> len(np,0);
  int end, begin;
  for(int ii=0; ii<np; ii++){
    end = double(N)/np*(ii+1);
    begin = double(N)/np*ii;
    len[ii] = 3*(end-begin);
  } 
  /*Position, Momemtum and Force arrays*/
  std::vector<double> Pos(len[id],0.0);
  std::vector<double> Mom(len[id],0.0);
  std::vector<double> Force(len[id],0.0);
  /*Fills position and momentum vectors with the initial conditions*/
  Initial_configuration(Pos, Mom, len, N,  tag, id, np, status);
  /*Caclculates total force felt by all particles*/
  Total_Force(Pos, Force, len, N, tag, id, np, status);
  /*Prints the force felt by all particles*/
  print_all(Force, len, N, tag, id, np, status);
  /*Finalizes MPI*/
  MPI_Finalize();
  return 0;
}

void Initial_configuration(std::vector<double>& Pos, std::vector<double>&Mom, std::vector<int>&len, int N, int tag, int id, int np, MPI_Status status){
  /*Process 0 creates the initial configuration of particles randomly
    (position and momentum) and then sends it to the other processes 
    according to their id*/
  if(id == 0){
    std::mt19937 gen(0);
    std::uniform_real_distribution<double> particles(0,1);
  
    std::vector<double> NPos(3*N);
    std::vector<double> NMom(3*N);
    for(int ii=0; ii<3*N;ii++){
      NPos[ii] = particles(gen);
      NMom[ii] = particles(gen);
      if(ii<len[0]){
	Pos[ii]=NPos[ii];
	Mom[ii]=NMom[ii];
      }
    }
    int from;
    for(int ii=1; ii<np; ii++){
      from = double(N)/np*(ii);
      MPI_Send(&NPos[3*from], len[ii], MPI_DOUBLE, ii , tag, MPI_COMM_WORLD);
      MPI_Send(&NMom[3*from], len[ii], MPI_DOUBLE, ii , tag, MPI_COMM_WORLD);  
    }
  }
  else{
    MPI_Recv(&Pos[0], len[id], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    MPI_Recv(&Mom[0], len[id], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
  }
}

void Total_Force(std::vector<double>& Pos, std::vector<double>& Force, std::vector<int>& len, int N, int tag, int id, int np, MPI_Status status){
  /*Temporal arrays for saving the position of particles shared in the ring*/
  int max = 3*N/np+3;
  std::vector<double> Temp(max,0.0);
  for(int ii=0; ii<len[id];ii++){
    Temp[ii] = Pos[ii];
  }
  std::vector<double> Temp2 = Temp;
  /*Each process calculates force due to its own particles*/
  Gravitational_force(Force,Pos,Temp,len[id]/3,len[id]/3);
  /*Ring*/ 
  int dst= (id+1)%np;
  int scr= (id-1+np)%np;
  for (int jj=0; jj<np-1; jj++){
    if (id%2==0){
      MPI_Send(&Temp[0], max, MPI_DOUBLE, dst , tag, MPI_COMM_WORLD);
      MPI_Recv(&Temp[0], max, MPI_DOUBLE, scr, tag, MPI_COMM_WORLD, &status);
      Gravitational_force(Force,Pos,Temp,len[id]/3,len[(scr-jj+np)%np]/3);
    }
    else{
      MPI_Recv(&Temp2[0], max, MPI_DOUBLE, scr, tag, MPI_COMM_WORLD, &status);
      MPI_Send(&Temp[0], max, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
      Gravitational_force(Force,Pos,Temp2,len[id]/3,len[(scr-jj+np)%np]/3);
      Temp = Temp2;
    }
  }
}

void Gravitational_force(std::vector<double>& Force, const std::vector<double>& vec0, const std::vector<double>& vec1, int len0, int len1){
  /*Calcultes force on particles in vec0 due to partcles in vec1*/
  double G, d2=0, pi=3.1415;
  G=4*pi*pi;
  for(int ii=0; ii<len0; ii++){
    for(int jj=0; jj<len1; jj++){
      d2=pow(vec0[3*ii+0]-vec1[3*jj+0],2)+pow(vec0[3*ii+1]-vec1[3*jj+1],2)+pow(vec0[3*ii+2]-vec1[3*jj+2],2);
      if (d2<1.0E-7){
	d2=1.0E-7;
      }
      for(int kk=0; kk<3; kk++){
	Force[3*ii+kk]+=G*(vec0[3*ii+kk]-vec1[3*jj+kk])/pow(d2,1.5);
      }
    }
  }
}

void print_all(std::vector<double>& Vec, const std::vector<int>& len, int N, int tag, int id, int np, MPI_Status status){
  /*Prints in terminal from process 0 the force felt by all particles*/
  if(id==0){
    std::cout.precision(2);
    std::cout<<std::scientific;
    print_vec(Vec, len[0]);
    std::vector<double> Temp(3*N/np+3,0.0);
    for (int kk =1; kk < np; kk++){
      MPI_Recv(&Temp[0], len[kk], MPI_DOUBLE, kk, tag, MPI_COMM_WORLD, &status);
      print_vec(Temp, len[kk]);   
    } 
  }
  else{
    MPI_Send(&Vec[0], len[id], MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
  }
}

void print_vec(std::vector<double>& Vec, int len){
  /*Prints each process force*/
  for (int ii = 0; ii < len; ii+=3){
    std::cout<<Vec[ii]<<"\t"<<Vec[ii+1]<<"\t"<<Vec[ii+2]<<"\n";
  }
}
