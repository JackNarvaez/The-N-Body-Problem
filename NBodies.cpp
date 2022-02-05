/*-----------------------------------------------
Evolution of a system of particles enclosed in a
box of unit sides, the only interaction is 
gravitation. The code implements the ring method.
---------------------------------------------*/

#include <mpi.h>
#include <iostream>
#include <fstream> 
#include <sstream> 
#include <vector>
#include <random>
#include <cstdlib>
#include <cmath>

void read_NParticles(const std::string &File_address, int& N){
  /*---------------------------------------------------------------------------
  read_data:
  Read number of particles from File_address.
  The data is supposed to be in the first line such as #  N.
  -----------------------------------------------------------------------------
  Arguments:
    File_address: File address from which the data is read.
    N : Number of particles.
  ---------------------------------------------------------------------------*/

  std::ifstream File;
  File.open (File_address, std::ifstream::in);    // Open file for reading
  std::string line;
  std::getline(File,line);
  std::istringstream iss(line);   // Separate line in columns
  std::string data;
  iss >> data;  // Delete sign #
  iss >> data;
  N = std::atoi(data.c_str());
  File.close();
}

void read_data(const std::string &File_address, std::vector<double>& Pos, std::vector<double>& Mom){
  /*---------------------------------------------------------------------------
  read_data:
  Read data about position, momentum and mass of particles.
  -----------------------------------------------------------------------------
  Arguments:
    File_address: File address from which the data is read.
    Pos   :   Position of particles (1D vector).
    Mom   :   Momentum of particles (1D vector).
  ---------------------------------------------------------------------------*/

  std::ifstream File;
  File.open (File_address, std::ifstream::in);    // Open file for reading
  std::string line;
	while (!File.eof()){
	  std::getline(File,line);
    // Omit empty lines and comments #
	  if (line.length() == 0 || line[0] == '#'){
	    continue;
    }else{
      std::istringstream iss(line);   // Separate line in columns
      std::string data;       
      for (int ii=0; ii < 3; ii++){
        iss >> data;
        Pos.push_back(atof(data.c_str()));  // Position
      }
      for (int ii=0; ii < 3; ii++){
        iss >> data;
        Mom.push_back(atof(data.c_str()));  // Momentum
      }
      iss >> data;
      Pos.push_back(atof(data.c_str()));  // Mass
    }
    }
    File.close();
}

void Initial_state(const std::string &File_address, std::vector<double>& Pos, std::vector<double>&Mom, std::vector<int>&len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
  /*---------------------------------------------------------------------------
  Initial_state:
  Process <root> reads the information about the initial configuration of the 
  system, fills the vectors <Pos> and <Mom>, and finally it distributes data
  among the other processes according to their <pId>.
  -----------------------------------------------------------------------------
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
    read_data(File_address, NPos, NMom);  // Read data
    // Save local information
    for(int ii=0; ii<4*len[root];ii++){
	    Pos[ii]=NPos[ii];
    }
    for(int ii=0; ii<3*len[root];ii++){
	    Mom[ii]=NMom[ii];
    }
    // Distribute data
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

void Gravitational_force(std::vector<double>& Force, const std::vector<double>& vec0, const std::vector<double>& vec1, const int & len0, const int & len1){
  /*---------------------------------------------------------------------------
  Gravitational_force:
  Calculate the gravitational force on particles in <vec0> due to particles in 
  <vec1>.
  -----------------------------------------------------------------------------
  Arguments:
    Force :   Force of particles in vec0 (1D vector)
    vec0  :   Local particles.
    vec1  :   Shared particles in the ring.
    len0  :   Number of local particles.
    len1  :   Size of vec1.
  ---------------------------------------------------------------------------*/

  const double G= 4*pow(M_PI,2);  // Gravitational constant [Msun*AU]
  double d2;  // Square of distance
  for(int ii=0; ii<len0; ii++){
    for(int jj=0; jj<len1; jj++){
      d2=pow(vec1[4*jj]-vec0[4*ii],2)+pow(vec1[4*jj+1]-vec0[4*ii+1],2)+pow(vec1[4*jj+2]-vec0[4*ii+2],2);
      if(d2<1.0E-10) d2=1.0E-7; // Lower distances are not valid
      for(int kk=0; kk<3; kk++){
	      Force[3*ii+kk]+=G*(vec1[4*jj+3])*(vec1[4*jj+kk]-vec0[4*ii+kk])/pow(d2,1.5);
      }
    }
  }
}

void Total_Force(const std::vector<double>& Pos, std::vector<double>& Force, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
  /*---------------------------------------------------------------------------
  Total_Force:
  Calculate the gravitational force for each particle due to all others.
  -----------------------------------------------------------------------------
  Arguments:
    Pos   :   Position of local particles (1D vector).
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

  //  Ring method 
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

void Save_vec(std::ofstream& File, const std::vector<double>& Vec, const int & len){
  /*---------------------------------------------------------------------------
  Save_vec:
  Save a vector <Vec> in Force.txt.
  -----------------------------------------------------------------------------
  Arguments:
    File  :   File where data is saved.
    Vec   :   Vector.
    N     :   Size of N.
  ---------------------------------------------------------------------------*/

  for (int ii = 0; ii < len; ii++){
    File << Vec[4*ii]<< "\t" << Vec[4*ii+1] << "\t" << Vec[4*ii+2] << std::endl;
  }
}

void Save_data(std::ofstream& File, const std::vector<double>& Pos, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
  /*---------------------------------------------------------------------------
  Save_Data:
  Save the Position of all particles in Evolution.txt.
  -----------------------------------------------------------------------------
  Arguments:
    File  :   File where data is saved.
    Pos   :   Position of particles (1D vector)
    len   :   Number of particles to send to each process.
    N     :   Number of particles.
    tag   :   Message tag.
    pId   :   Process identity.
    nP    :   Number of processes.
    root  :   Root process which reads data.
    status:   Status object.
  ---------------------------------------------------------------------------*/

  // Collect results in <root> process.
  if(pId==root){
    std::cout.precision(2);
    std::cout<<std::scientific;
    File.open("Evolution.txt",std::fstream::app);
    Save_vec(File, Pos, len[root]);
    std::vector<double> Temp(4*(N/nP+1),0.0);
    for (int ii =0; ii < nP; ii++){
      if (ii != pId){
        MPI_Recv(&Temp[0], 4*len[ii], MPI_DOUBLE, ii, tag, MPI_COMM_WORLD, &status);
        Save_vec(File, Temp,  len[ii]);
      }  
    }
    File.close(); 
  }else{
    MPI_Send(&Pos[0], 4*len[pId], MPI_DOUBLE, root, tag, MPI_COMM_WORLD);
  }
}

void Euler(std::vector<double>& Pos, std::vector<double>& Mom, const std::vector<double>& Force, const double &dt, const int & N){
  /*---------------------------------------------------------------------------
  Euler:
  Euler method to calculate next position and momentum.
  -----------------------------------------------------------------------------
  Arguments:

    Pos   :   Position of particles (1D vector)
    Mom   :   Momentum of particles (1D vector)
    Force :   Force of particles in vec0 (1D vector)
    N     :   Local particles for each process.
  ---------------------------------------------------------------------------*/
  for (int ii=0; ii < N; ii++){
    for (int jj=0; jj<3; jj++){
      Mom[3*ii+jj] += Force[3*ii+jj]*dt;
      Pos[4*ii+jj] += Mom[3*ii+jj]*dt;
    }
  }
}

void RK45(std::vector<double>& Pos, std::vector<double>& Mom, const std::vector<double>& Force, const double &dt, const int & N){
  /*---------------------------------------------------------------------------
  Euler:
  Euler method to calculate next position and momentum.
  -----------------------------------------------------------------------------
  Arguments:

    Pos   :   Position of particles (1D vector)
    Mom   :   Momentum of particles (1D vector)
    Force :   Force of particles in vec0 (1D vector)
    N     :   Local particles for each process.
  ---------------------------------------------------------------------------*/
  for (int ii=0; ii < N; ii++){
    for (int jj=0; jj<3; jj++){
      Mom[3*ii+jj] += Force[3*ii+jj]*dt;
      Pos[4*ii+jj] += Mom[3*ii+jj]/Pos[4*ii+3]*dt;
    }
  }
}

void Evolution(std::ofstream& File, std::vector<double>& Pos, std::vector<double>& Mom, std::vector<double>& Force, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status, const int &steps, const double & dt, const int & jump){
  /*---------------------------------------------------------------------------
  Evolution:
  Evolution of the system of particles under gravitational interactions.
  -----------------------------------------------------------------------------
  Arguments:
    File  :   File where data is saved.
    Pos   :   Position of particles (1D vector)
    Mom   :   Momentum of particles (1D vector)
    Force :   Force of particles in vec0 (1D vector)
    len   :   Number of particles to send to each process.
    N     :   Number of particles.
    tag   :   Message tag.
    pId   :   Process identity.
    nP    :   Number of processes.
    root  :   Root process which reads data.
    status:   Status object.
    steps :   Number of steps
    dt    :   Size of step
    jump  :   Number of jumped steps
  ---------------------------------------------------------------------------*/
  
  if (pId==root){
    remove("Evolution.txt");
  }
  Save_data(File, Pos, len, N, tag, pId, nP, root, status);
  for (int ii=0; ii<steps; ii++){
    Total_Force(Pos, Force, len, N, tag, pId, nP, root, status);
    Euler(Pos, Mom, Force, dt, len[pId]);
    if ( ii%jump == 0) Save_data(File, Pos, len, N, tag, pId, nP, root, status);  
  }
}