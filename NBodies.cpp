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

using function = void(const std::vector<double>&, std::vector<double>&, const std::vector<int>&, const int &, const int &, const int &, const int &, const int &, MPI_Status);

void read_data(const std::string &File_address, std::vector<double>& Pos, std::vector<double>& Vel){
  /*---------------------------------------------------------------------------
  read_data:
  Read data about position, velocity and mass of particles.
  -----------------------------------------------------------------------------
  Arguments:
    File_address: File address from which the data is read.
    Pos   :   Position of particles (1D vector).
    Vel   :   Velocity of particles (1D vector).
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
      Pos.push_back(atof(data.c_str()));  // Mass
    }
    }
    File.close();
}

void Gravitational_Acc(std::vector<double>& Acc, const std::vector<double>& vec0, const std::vector<double>& vec1, const int & len0, const int & len1){
  /*---------------------------------------------------------------------------
  Gravitational_Acc:
  Calculate the gravitational acceleration on particles in <vec0> due to particles in 
  <vec1>.
  -----------------------------------------------------------------------------
  Arguments:
    Acc   :   Acceleration of particles in vec0 (1D vector)
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
	      Acc[3*ii+kk]+=G*(vec1[4*jj+3])*(vec1[4*jj+kk]-vec0[4*ii+kk])/pow(d2,1.5);
      }
    }
  }
}

void Acceleration(const std::vector<double>& Pos, std::vector<double>& Acc, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
  /*---------------------------------------------------------------------------
  Acceleration:
  Calculate the gravitational acceleration for each particle due to all others.
  -----------------------------------------------------------------------------
  Arguments:
    Pos   :   Position of local particles (1D vector).
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
  int max = 4*(N/nP+1);
  std::vector<double> Temp(max,0.0);
  for(int ii=0; ii<4*len[pId]; ii++){
    Temp[ii] = Pos[ii];
  }
  std::vector<double> Temp2(Temp);
  
  //  Gravitational acceleration due to local particles
  Gravitational_Acc(Acc, Pos, Temp, len[pId], len[pId]);

  //  Ring method 
  int dst= (pId+1)%nP;
  int scr= (pId-1+nP)%nP;
  for (int jj=0; jj<nP-1; jj++){
    if (pId%2==0){
      MPI_Send(&Temp[0], max, MPI_DOUBLE, dst , tag, MPI_COMM_WORLD);
      MPI_Recv(&Temp[0], max, MPI_DOUBLE, scr, tag, MPI_COMM_WORLD, &status);
      Gravitational_Acc(Acc, Pos, Temp, len[pId], len[(scr-jj+nP)%nP]);
    }
    else{
      MPI_Recv(&Temp2[0], max, MPI_DOUBLE, scr, tag, MPI_COMM_WORLD, &status);
      MPI_Send(&Temp[0], max, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
      Gravitational_Acc(Acc, Pos, Temp2, len[pId], len[(scr-jj+nP)%nP]);
      Temp = Temp2;
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

void Euler(std::vector<double>& Pos, std::vector<double>& Vel, std::vector<double>& Acc, const double &dt, const int & N, function Accel, const std::vector<int>& len, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
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
  Accel(Pos, Acc, len, N, tag, pId, nP, root, status);
  for (int ii=0; ii < len[pId]; ii++){
    for (int jj=0; jj<3; jj++){
      Vel[3*ii+jj] += Acc[3*ii+jj]*dt;
      Pos[4*ii+jj] += Vel[3*ii+jj]*dt;
    }
  }
}

void PEFRL(std::vector<double>& Pos, std::vector<double>& Vel, std::vector<double>& Acc, const double &dt, const int & N, function Accel, const std::vector<int>& len, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
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
      X[4*ii+jj] += xi*dt*V[3*ii+jj];
    }
  }
  std::fill (Acc.begin(),Acc.end(),0);
  Accel(X, Acc, len, N, tag, pId, nP, root, status);
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      V[3*ii+jj] +=  0.5*(1.-2*gamma)*dt*Acc[3*ii+jj];
    }
  }
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      X[4*ii+jj] += chi*dt*V[3*ii+jj];
    }
  }
  std::fill (Acc.begin(),Acc.end(),0);
  Accel(X, Acc, len, N, tag, pId, nP, root, status);
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      V[3*ii+jj] += gamma*dt*Acc[3*ii+jj];
    }
  }
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      X[4*ii+jj] += (1.-2*(chi+xi))*dt*V[3*ii+jj];
    }
  }
  std::fill (Acc.begin(),Acc.end(),0);
  Accel(X, Acc, len, N, tag, pId, nP, root, status);
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      V[3*ii+jj] += gamma*dt*Acc[3*ii+jj];
    }
  }
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      X[4*ii+jj] += chi*dt*V[3*ii+jj];
    }
  }
  std::fill (Acc.begin(),Acc.end(),0);
  Accel(X, Acc, len, N, tag, pId, nP, root, status);
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      Vel[3*ii+jj] = V[3*ii+jj]+0.5*(1.-2*gamma)*dt*Acc[3*ii+jj];
    }
  }
  for (int ii = 0; ii < local_particles; ii++){
    for (int jj = 0; jj < 3; jj++){
      Pos[4*ii+jj] = X[4*ii+jj]+xi*dt*Vel[3*ii+jj];
    }
  }
}

void Evolution(std::ofstream& File, std::vector<double>& Pos, std::vector<double>& Vel, std::vector<double>& Acc, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status, const int &steps, const double & dt, const int & jump, function Accel){
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
    //Euler(Pos, Vel, Acc, dt, N, Accel, len, tag, pId, nP, root, status);
    PEFRL(Pos, Vel, Acc, dt, N, Accel, len, tag, pId, nP, root, status);
    if ( ii%jump == 0) Save_data(File, Pos, len, N, tag, pId, nP, root, status);
    std::fill (Acc.begin(),Acc.end(),0);
  }
}