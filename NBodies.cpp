/*-----------------------------------------------
Evolution of a system of particles. The only 
interaction is gravitation. 
The code implements the ring method.
-----------------------------------------------*/

#include <mpi.h>
#include <iostream>
#include <fstream> 
#include <sstream> 
#include <vector>
#include <random>
#include <cstdlib>
#include <cmath>
#include <string>

using function = void(const std::vector<double>&, const std::vector<double>&, std::vector<double>&, const int *, const int &, const int &, const int &, const int &, const int &, MPI_Status);
using Integrator = void(std::vector<double>&, std::vector<double>&, const std::vector<double>&, std::vector<double>&, const double &, const int &, function, const int *, const int &, const int &, const int &, const int &, MPI_Status);

void read_data(const std::string &File_address, std::vector<double>& Pos, std::vector<double>& Vel, std::vector<double>& Mass){
  /*---------------------------------------------------------------------------
  read_data:
  Read data about position, velocity and mass of particles.
  -----------------------------------------------------------------------------
  Arguments:
    File_address: File address from which the data is read.
    Pos   :   Position of particles (1D vector).
    Vel   :   Velocity of particles (1D vector).
    Mass  :   Mass of particles (1D vector).
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
        Vel.push_back(atof(data.c_str()));  // Velocity
      }
      iss >> data;
      Mass.push_back(atof(data.c_str()));  // Mass
    }
    }
    File.close();
}

void Gravitational_Acc(std::vector<double>& Acc, const std::vector<double>& Pos0, const std::vector<double>& Pos1, const std::vector<double>& Mass1, const int & len0, const int & len1){
  /*---------------------------------------------------------------------------
  Gravitational_Acc:
  Calculate the gravitational acceleration on particles in <Pos0> due to particles in 
  <Pos1>.
  -----------------------------------------------------------------------------
  Arguments:
    Acc   :   Acceleration of particles in vec0 (1D vector)
    Pos0  :   Local particles.
    Pos1  :   Shared particles in the ring.
    Mass1 :   Mass of shared particles.
    len0  :   Number of local particles.
    len1  :   Size of vec1.
  ---------------------------------------------------------------------------*/

  const double G= 4*pow(M_PI,2);  // Gravitational constant [Msun*AU]
  double d2;  // Square of distance
  for(int ii=0; ii<len0; ii++){
    for(int jj=0; jj<len1; jj++){
      d2=pow(Pos1[3*jj]-Pos0[3*ii],2)+pow(Pos1[3*jj+1]-Pos0[3*ii+1],2)+pow(Pos1[3*jj+2]-Pos0[3*ii+2],2);
      if(d2<1.0E-10) d2=1.0E-7; // Lower distances are not valid
      for(int kk=0; kk<3; kk++){
	      Acc[3*ii+kk]+=G*(Mass1[jj])*(Pos1[3*jj+kk]-Pos0[3*ii+kk])/pow(d2,1.5);
      }
    }
  }
}

void Acceleration(const std::vector<double>& Pos, const std::vector<double>& Mass, std::vector<double>& Acc, const int * len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
  /*---------------------------------------------------------------------------
  Acceleration:
  Calculate the gravitational acceleration for each particle due to all others.
  -----------------------------------------------------------------------------
  Arguments:
    Pos   :   Position of local particles (1D vector).
    Mass  :   Mass of local particles (1D vector)
    Acc   :   Acceleration of local particles (1D vector)
    len   :   Number of particles to send to each process.
    N     :   Number of total particles.
    tag   :   Message tag.
    pId   :   Process identity.
    nP    :   Number of processes.
    root  :   Root process which reads data.
    status:   Status object.
  ---------------------------------------------------------------------------*/

  //  Temporal arrays for saving data of shared particles along the ring
  int max = (N/nP+1);
  std::vector<double> BufferPos(Pos);
  std::vector<double> BufferMass(Mass);
  BufferPos.insert(BufferPos.end(), 3, 0);
  BufferMass.push_back(0);
  std::vector<double> Buffer2Pos(3*max, 0);
  std::vector<double> Buffer2Mass(max, 0);
  
  //  Gravitational acceleration due to local particles
  Gravitational_Acc(Acc, Pos, BufferPos, Mass, len[pId], len[pId]);
  
  //  Ring method 
  int dst= (pId+1)%nP;
  int scr= (pId-1+nP)%nP;
  for (int jj=0; jj<nP-1; jj++){
    if (pId%2==0){
      MPI_Send(&BufferPos[0], 3*max, MPI_DOUBLE, dst , tag, MPI_COMM_WORLD);
      MPI_Recv(&BufferPos[0], 3*max, MPI_DOUBLE, scr, tag, MPI_COMM_WORLD, &status);
      MPI_Send(&BufferMass[0], max, MPI_DOUBLE, dst , tag, MPI_COMM_WORLD);
      MPI_Recv(&BufferMass[0], max, MPI_DOUBLE, scr, tag, MPI_COMM_WORLD, &status);
      Gravitational_Acc(Acc, Pos, BufferPos, BufferMass, len[pId], len[(scr-jj+nP)%nP]);
    }
    else{
      MPI_Recv(&Buffer2Pos[0], 3*max, MPI_DOUBLE, scr, tag, MPI_COMM_WORLD, &status);
      MPI_Send(&BufferPos[0], 3*max, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
      MPI_Recv(&Buffer2Mass[0], max, MPI_DOUBLE, scr, tag, MPI_COMM_WORLD, &status);
      MPI_Send(&BufferMass[0], max, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
      Gravitational_Acc(Acc, Pos, Buffer2Pos, Buffer2Mass, len[pId], len[(scr-jj+nP)%nP]);
      BufferPos = Buffer2Pos;
      BufferMass = Buffer2Mass;
    }
  }
}

void Save_vec(std::ofstream& File, const std::vector<double>& Vec, const int & N){
  /*---------------------------------------------------------------------------
  Save_vec:
  Save a vector <Vec> in File.
  -----------------------------------------------------------------------------
  Arguments:
    File  :   File where data is saved.
    Vec   :   Vector.
    N     :   Size of N.
  ---------------------------------------------------------------------------*/

  for (int ii = 0; ii < N; ii++){
    File << Vec[3*ii]<< "\t" << Vec[3*ii+1] << "\t" << Vec[3*ii+2] << std::endl;
  }
}

void Save_data(std::ofstream& File, const std::vector<double>& Pos, const int * len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
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
    std::vector<double> Temp(3*(N/nP+1),0.0);
    for (int ii =0; ii < nP; ii++){
      if (ii != pId){
        MPI_Recv(&Temp[0], 3*len[ii], MPI_DOUBLE, ii, tag, MPI_COMM_WORLD, &status);
        Save_vec(File, Temp,  len[ii]);
      }  
    }
    File.close(); 
  }else{
    MPI_Send(&Pos[0], 3*len[pId], MPI_DOUBLE, root, tag, MPI_COMM_WORLD);
    }
}

void Euler(std::vector<double>& Pos, std::vector<double>& Vel, const std::vector<double>& Mass, std::vector<double>& Acc, const double &dt, const int & N, function Accel, const int * len, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
  /*---------------------------------------------------------------------------
  Euler:
  Euler method to calculate next position and velocity.
  -----------------------------------------------------------------------------
  Arguments:
    Pos   :   Position of particles (1D vector)
    Vel   :   Velocity of particles (1D vector)
    Acc   :   Acceleration of particles in vec0 (1D vector)
    N     :   Local particles for each process.
  ---------------------------------------------------------------------------*/
  Accel(Pos, Mass, Acc, len, N, tag, pId, nP, root, status);
  for (int ii=0; ii < len[pId]; ii++){
    for (int jj=0; jj<3; jj++){
      Vel[3*ii+jj] += Acc[3*ii+jj]*dt;
      Pos[3*ii+jj] += Vel[3*ii+jj]*dt;
    }
  }
}

void PEFRL(std::vector<double>& Pos, std::vector<double>& Vel, const std::vector<double>& Mass, std::vector<double>& Acc, const double &dt, const int & N, function Accel, const int * len, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
  /*---------------------------------------------------------------------------
  Euler:
  Euler method to calculate next position and velocity.
  -----------------------------------------------------------------------------
  Arguments:
    Pos   :   Position of particles (1D vector)
    Vel   :   Velocity of particles (1D vector)
    Acc   :   Acceleration of particles in vec0 (1D vector)
    N     :   Local particles for each process.
  ---------------------------------------------------------------------------*/
  int local_particles = len[pId];
  std::vector<double> X(Pos);
  std::vector<double> V(Vel);

  // Parameters
  double xi = 0.1786178958448091e+0;
  double gamma = -0.2123418310626054e+0;
  double chi = -0.6626458266981849e-1;

  // Main loop
  for(int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      X[3*ii+jj] += xi*dt*V[3*ii+jj];
    }
  }
  std::fill (Acc.begin(),Acc.end(),0);
  Accel(X, Mass, Acc, len, N, tag, pId, nP, root, status);
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      V[3*ii+jj] +=  0.5*(1.-2*gamma)*dt*Acc[3*ii+jj];
    }
  }
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      X[3*ii+jj] += chi*dt*V[3*ii+jj];
    }
  }
  std::fill (Acc.begin(),Acc.end(),0);
  Accel(X, Mass, Acc, len, N, tag, pId, nP, root, status);
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      V[3*ii+jj] += gamma*dt*Acc[3*ii+jj];
    }
  }
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      X[3*ii+jj] += (1.-2*(chi+xi))*dt*V[3*ii+jj];
    }
  }
  std::fill (Acc.begin(),Acc.end(),0);
  Accel(X, Mass, Acc, len, N, tag, pId, nP, root, status);
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      V[3*ii+jj] += gamma*dt*Acc[3*ii+jj];
    }
  }
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      X[3*ii+jj] += chi*dt*V[3*ii+jj];
    }
  }
  std::fill (Acc.begin(),Acc.end(),0);
  Accel(X, Mass, Acc, len, N, tag, pId, nP, root, status);
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      Vel[3*ii+jj] = V[3*ii+jj]+0.5*(1.-2*gamma)*dt*Acc[3*ii+jj];
    }
  }
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      Pos[3*ii+jj] = X[3*ii+jj]+xi*dt*Vel[3*ii+jj];
    }
  }
}

void Evolution(std::ofstream& File, std::vector<double>& Pos, std::vector<double>& Vel, const std::vector<double>& Mass, std::vector<double>& Acc, const int * len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status, const int &steps, const double & dt, const int & jump, function Accel, Integrator evol){
  /*---------------------------------------------------------------------------
  Evolution:
  Evolution of the system of particles under gravitational interactions.
  -----------------------------------------------------------------------------
  Arguments:
    File  :   File where data is saved.
    Pos   :   Position of particles (1D vector)
    Vel   :   Velocity of particles (1D vector)
    Acc   :   Acceleration of particles in vec0 (1D vector)
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
    evol(Pos, Vel, Mass, Acc, dt, N, Accel, len, tag, pId, nP, root, status);
    if ( ii%jump == 0) Save_data(File, Pos, len, N, tag, pId, nP, root, status);
    std::fill (Acc.begin(),Acc.end(),0);
  }
}