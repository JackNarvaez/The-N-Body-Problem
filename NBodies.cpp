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

void read_data(const std::string &File_address, std::vector<double>& Pos, std::vector<double>& Mom){
  /*-------------------------------------------------------------------------
  read_data:
  Read data about position, momentum and mass of particles.
  ---------------------------------------------------------------------------
  Arguments:
    File_address: File address from which the data is read.
    Pos   :   Position of particles (1D vector).
    Mom   :   Momentum of particles (1D vector).
  ---------------------------------------------------------------------------*/
  std::ifstream File;
  File.open (File_address, std::ifstream::in);    // Open file
  std::string line;
	while (!File.eof()){
	  std::getline(File,line);
    // Omit empty lines and comments
	  if (line.length() == 0 || line[0] == '#'){
	    continue;
    }else{
      std::istringstream iss(line);   // Separate line in columns
      std::string data;       
      for (int ii=0; ii < 3; ii++){
        iss >> data;
        Pos.push_back(atof(data.c_str()));
      }
      for (int ii=0; ii < 3; ii++){
        iss >> data;
        Mom.push_back(atof(data.c_str()));
      }
      iss >> data;
      Pos.push_back(atof(data.c_str())); 
    }
    }
    File.close();
}

void Initial_state(const std::string &File_address, std::vector<double>& Pos, std::vector<double>&Mom, std::vector<int>&len, int N, int tag, int pId, int nP, int root, MPI_Status status){
  /*-------------------------------------------------------------------------
  Initial_state:
  Process root reads the information about the initial configuration of the 
  system, fills the data Pos, Mom and Mass, and finally it distribuites data
  among the other processes according to their id
  ---------------------------------------------------------------------------
  Arguments:
    File_address: File address from which the data is read.
    Pos   :   Position and mass of local particles (1D vector).
    Mom   :   Momentum of local particles (1D vector).
    len   :   Number of particles to send to each process.
    N     :   Number of total particles.
    tag   :   Message tag.
	  pId   :   Process identity.
    nP    :   Number of processes.
    root  :   Root process which reads data.
    status:   Status object.
  ---------------------------------------------------------------------------*/
  
  if(pId == root){
    std::vector<double> NPos;
    std::vector<double> NMom;
    read_data(File_address, NPos, NMom);
    for(int ii=0; ii<4*len[root];ii++){
	    Pos[ii]=NPos[ii];
    }
    for(int ii=0; ii<3*len[root];ii++){
	    Mom[ii]=NMom[ii];
    }
    int Lower;
    for(int ii=0; ii<nP; ii++){
      if (ii != root){
        Lower = double(N)/nP*(ii);
        MPI_Send(&NPos[4*Lower], 4*len[ii], MPI_DOUBLE, ii , tag, MPI_COMM_WORLD);
        MPI_Send(&NMom[3*Lower], 3*len[ii], MPI_DOUBLE, ii , tag, MPI_COMM_WORLD);
      }
    }
  }else{
    MPI_Recv(&Pos[0], 4*len[pId], MPI_DOUBLE, root, tag, MPI_COMM_WORLD, &status);
    MPI_Recv(&Mom[0], 3*len[pId], MPI_DOUBLE, root, tag, MPI_COMM_WORLD, &status);
  }
}

void Gravitational_force(std::vector<double>& Force, const std::vector<double>& vec0, const std::vector<double>& vec1, int len0, int len1){
  /*-------------------------------------------------------------------------
  Gravitational_force:
  Calculate the gravitational force on particles in vec0 due to particles in 
  vec1.
  ---------------------------------------------------------------------------
  Arguments:
    Force :   Force of particles in vec0 (1D vector)
    vec0  :   Local particles.
    vec1  :   Shared particles in the ring.
    len0  :   Number of local particles.
    len1  :   Size of vec1.
  ---------------------------------------------------------------------------*/
  const double G= 4*pow(M_PI,2);
  double d2;
  for(int ii=0; ii<len0; ii++){
    for(int jj=0; jj<len1; jj++){
      d2=pow(vec1[4*jj]-vec0[4*ii],2)+pow(vec1[4*jj+1]-vec0[4*ii+1],2)+pow(vec1[4*jj+2]-vec0[4*ii+2],2);
      if(d2<1.0E-10)d2=1.0E-7;
      for(int kk=0; kk<3; kk++){
	      Force[3*ii+kk]+=G*(vec0[4*ii+3]*vec1[4*ii+3])*(vec1[4*jj+kk]-vec0[4*ii+kk])/pow(d2,1.5);
      }
    }
  }
}

void Total_Force(std::vector<double>& Pos, std::vector<double>& Force, std::vector<int>& len, int N, int tag, int pId, int nP, int root, MPI_Status status){
  /*-------------------------------------------------------------------------
  Total_Force:
  Calculate the gravitational force for each partcicle due to all others.
  ---------------------------------------------------------------------------
  Arguments:
    Pos   :   Position of local particles (1D vector).
    Mass  :   Mass of local particles (1D vector).
    Force :   Force of local particles (1D vector)
    len   :   Number of particles to send to each process.
    N     :   Number of total particles.
    tag   :   Message tag.
	  pId   :   Process identity.
    nP    :   Number of processes.
    root  :   Root process which reads data.
    status:   Status object.
  ---------------------------------------------------------------------------*/

  //  Temporal arrays for saving data of shared particles along the ring
  int max = 4*(N/nP+1);
  std::vector<double> Temp(max,0.0);
  for(int ii=0; ii<4*len[pId]; ii++){
    Temp[ii] = Pos[ii];
  }
  std::vector<double> Temp2(Temp);
  
  //  Gravitational force due to local particles
  Gravitational_force(Force, Pos, Temp, len[pId], len[pId]);

  //  Ring 
  int dst= (pId+1)%nP;
  int scr= (pId-1+nP)%nP;
  for (int jj=0; jj<nP-1; jj++){
    if (pId%2==0){
      MPI_Send(&Temp[0], max, MPI_DOUBLE, dst , tag, MPI_COMM_WORLD);
      MPI_Recv(&Temp[0], max, MPI_DOUBLE, scr, tag, MPI_COMM_WORLD, &status);
      Gravitational_force(Force, Pos, Temp, len[pId], len[(scr-jj+nP)%nP]);
    }
    else{
      MPI_Recv(&Temp2[0], max, MPI_DOUBLE, scr, tag, MPI_COMM_WORLD, &status);
      MPI_Send(&Temp[0], max, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
      Gravitational_force(Force, Pos, Temp2, len[pId], len[(scr-jj+nP)%nP]);
      Temp = Temp2;
    }
  }
}

void Save_vec(std::ofstream& File, std::vector<double>& Vec, int len){
  /*Save_vec:
  Save a vector Vec in Force.txt.
  ---------------------------------------------------------------------------
  Arguments:
    File  :   File where data is saved.
    Vec   :   Vector.
    N     :   Size of N.
  ---------------------------------------------------------------------------*/
  File.open("Force.txt", std::fstream::app);
  for (int ii = 0; ii < len; ii++){
    File << Vec[3*ii]<< "\t" << Vec[3*ii+1] << "\t" << Vec[3*ii+2] << std::endl;
  }
  File.close();
}

void Save_data(std::ofstream& File, std::vector<double>& Force, const std::vector<int>& len, int N, int tag, int pId, int nP, int root, MPI_Status status){
  /*-------------------------------------------------------------------------
  Save_Data:
  Save the gravitation forces of all particles in Force.txt.
  ---------------------------------------------------------------------------
  Arguments:
    File  :   File where data is saved.
    Force :   Force of particles (1D vector)
    len   :   Number of particles to send to each process.
    N     :   Number of particles.
    tag   :   Message tag.
    pId   :   Process identity.
    nP    :   Number of processes.
    root  :   Root process which reads data.
    status:   Status object.
  ---------------------------------------------------------------------------*/
  if(pId==root){
    std::cout.precision(2);
    std::cout<<std::scientific;
    Save_vec(File, Force, len[root]);
    std::vector<double> Temp(3*(N/nP+1),0.0);
    for (int ii =0; ii < nP; ii++){
      if (ii != pId){
        MPI_Recv(&Temp[0], 3*len[ii], MPI_DOUBLE, ii, tag, MPI_COMM_WORLD, &status);
        Save_vec(File, Temp,  len[ii]);
      }  
    } 
  }else{
    MPI_Send(&Force[0], 3*len[pId], MPI_DOUBLE, root, tag, MPI_COMM_WORLD);
  }
}