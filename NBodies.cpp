/*-----------------------------------------------------------------------------
The N-Body problem using MPI.
The code implements the ring method.
-----------------------------------------------------------------------------*/

#include <mpi.h>
#include <iostream>
#include <fstream> 
#include <sstream> 
#include <cstdlib>
#include <cmath>
#include <string>

using function = void(const double *, const double *, double *, const int *, const int &, const int &, const int &, const int &, const int &, MPI_Status);
using Integrator = void(double *, double *, const double *, double *, const double &, const int &, function, const int *, const int &, const int &, const int &, const int &, MPI_Status);

void read_data(const std::string & File_address, double * Pos, double * Vel, double * Mass) {
  /*---------------------------------------------------------------------------
  Reads data from bodies: position, velocity and mass.
  -----------------------------------------------------------------------------
  Arguments:
    File_address:   File address from where the data is read.
    Pos         :   Position (1D vector) [xi, yi, zi].
    Vel         :   Velocity (1D vector) [vxi, vyi, vzi].
    Mass        :   Mass (1D vector) [mi].
  ---------------------------------------------------------------------------*/

  std::ifstream File;
  File.open (File_address, std::ifstream::in);    // Open file for reading
  std::string line;
  int row = 0;
  int ii;
	while (!File.eof()){
	  std::getline(File,line);
    // Omit empty lines and comments #
	  if (line.length() == 0 || line[0] == '#'){
	    continue;
    }else{
      std::istringstream iss(line);               // Separate line in columns
      std::string data;       
      for (ii = 0; ii < 3; ii++){
        iss >> data;
        Pos[ii + 3*row] = atof(data.c_str());     // Position
      }
      for (ii = 0; ii < 3; ii++){
        iss >> data;
        Vel[ii + 3*row] = atof(data.c_str());     // Velocity
      }
      iss >> data;
      Mass[row] = atof(data.c_str());             // Mass
      row += 1;
    }
  }
  File.close();
}

void Gravitational_Acc(double * Acc, const double * Pos0, const double * Pos1, const double * Mass1, const int & len0, const int & len1){
  /*---------------------------------------------------------------------------
  Calculates the gravitational acceleration on bodies in Pos0 due to bodies in 
  Pos1.
  -----------------------------------------------------------------------------
  Arguments:
    Acc   :   Acceleration of bodies in Pos0 (1D vector).
    Pos0  :   Position of local bodies.
    Pos1  :   Position of shared bodies in the ring.
    Mass1 :   Mass of shared bodies.
    len0  :   Number of local bodies.
    len1  :   Size of Pos1.
  ---------------------------------------------------------------------------*/

  const double G= 4*pow(M_PI,2);                  // Gravitational constant [Msun*AU]
  int ii, jj;
  double drelx, drely, drelz, d2, inv_rtd2, cb_d2;
  for(ii = 0; ii < len0; ii++){
    for(jj = 0; jj < len1; jj++){
      drelx = Pos1[3*jj]-Pos0[3*ii];
      drely = Pos1[3*jj+1]-Pos0[3*ii+1];
      drelz = Pos1[3*jj+2]-Pos0[3*ii+2];
      d2 = drelx*drelx+drely*drely+drelz*drelz;
      if(d2<1.0E-10) d2=1.0E-7;                   // Lower distances are not valid
      inv_rtd2 = 1./sqrt(d2);
      cb_d2 = inv_rtd2*inv_rtd2*inv_rtd2;

	    Acc[3*ii]+=G*Mass1[jj]*drelx*cb_d2;
	    Acc[3*ii+1]+=G*Mass1[jj]*drely*cb_d2;
	    Acc[3*ii+2]+=G*Mass1[jj]*drelz*cb_d2;
    }
  }
}

void Acceleration(const double * Pos, const double * Mass, double * Acc, const int * len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
  /*---------------------------------------------------------------------------
  Calculates the gravitational acceleration for each body due to all others.
  -----------------------------------------------------------------------------
  Arguments:
    Pos   :   Position of local bodies (1D vector).
    Mass  :   Mass of local bodies (1D vector).
    Acc   :   Acceleration of local bodies (1D vector).
    len   :   Array with the number of bodies per node.
    N     :   Total bodies.
    tag   :   Message tag.
    pId   :   Process identity.
    nP    :   Number of processes.
    root  :   Root process.
    status:   Status object.
  ---------------------------------------------------------------------------*/

  int max = N/nP+1;

  //  Temporal arrays for saving data of shared bodies along the ring
  double *BufferPos = (double *) malloc(3*(len[pId]+1)*sizeof(double));
  double *BufferMass = (double *) malloc((len[pId]+1)*sizeof(double));

  for (int ii=0; ii < len[pId]; ii++){
    BufferPos[3*ii]   = Pos[3*ii];
    BufferPos[3*ii+1] = Pos[3*ii+1];
    BufferPos[3*ii+2] = Pos[3*ii+2];
    BufferMass[ii]    = Mass[ii];
  }

  BufferPos[3*len[pId]]   = 0.0;
  BufferPos[3*len[pId]+1] = 0.0;
  BufferPos[3*len[pId]+2] = 0.0;
  BufferMass[len[pId]]    = 0.0;
  
  //  Gravitational acceleration due to local particles
  Gravitational_Acc(Acc, Pos, Pos, Mass, len[pId], len[pId]);
  
  /*The Ring Loop*/
  int dst= (pId+1)%nP;
  int scr= (pId-1+nP)%nP;
  
  for (int jj=0; jj<nP-1; jj++){
    MPI_Sendrecv_replace(BufferPos, 3*max, MPI_DOUBLE, dst, tag, scr, tag, MPI_COMM_WORLD, &status);
    MPI_Sendrecv_replace(BufferMass, max, MPI_DOUBLE, dst, tag, scr, tag, MPI_COMM_WORLD, &status);
    Gravitational_Acc(Acc, Pos, BufferPos, BufferMass, len[pId], len[(scr-jj+nP)%nP]);
  }

  free(BufferPos);
  free(BufferMass);
}

void Save_vec(std::ofstream & File, const double * Pos, const double * Mass, const int & N){
  /*---------------------------------------------------------------------------
  Saves info from Pos and Mass in File.
  -----------------------------------------------------------------------------
  Arguments:
    File  :   File where data is saved.
    Pos   :   Positions.
    Mass  :   Masses.
    N     :   Size of Mass.
  ---------------------------------------------------------------------------*/
  int ii;
  for (ii = 0; ii < N; ii++){
    File << Pos[3*ii]<< "\t" << Pos[3*ii+1] << "\t" << Pos[3*ii+2]<< "\t" << Mass[ii] << std::endl;
  }
}

void Save_data(std::ofstream & File, const double * Pos, const double * Mass, const int * len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
  /*---------------------------------------------------------------------------
  Saves positions and masses of all bodies in Evolution.txt.
  -----------------------------------------------------------------------------
  Arguments:
    File  :   File where data is saved.
    Pos   :   Position of bodies (1D vector).
    Mass  :   Mass of bodies (1D vector).
    len   :   Array with the number of bodies per node.
    N     :   Number of bodies.
    tag   :   Message tag.
    pId   :   Process identity.
    nP    :   Number of processes.
    root  :   Root process.
    status:   Status object.
  ---------------------------------------------------------------------------*/

  // Collect results in <root> process.
  if(pId==root){
    std::cout.precision(2);
    std::cout<<std::scientific;
    File.open("Evolution.txt",std::fstream::app);
    Save_vec(File, Pos, Mass, len[root]);
    double * TempP = (double *) malloc(3*(N/nP+1)*sizeof(double));
    double * TempM = (double *) malloc((N/nP+1)*sizeof(double));
    int ii;
    for (ii = 0; ii < N/nP+1; ii++) {
      TempP[3*ii] = 0.0;
      TempP[3*ii+1] = 0.0;
      TempP[3*ii+2] = 0.0;
      TempM[ii] = 0.0;
    }

    for (ii = 0; ii < nP; ii++){
      if (ii != pId){
        MPI_Recv(TempP, 3*len[ii], MPI_DOUBLE, ii, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(TempM, len[ii], MPI_DOUBLE, ii, tag, MPI_COMM_WORLD, &status);
        Save_vec(File, TempP, TempM,  len[ii]);
      }
    }
    File.close();
    free(TempP);
    free(TempM);
  }else{
    MPI_Send(Pos, 3*len[pId], MPI_DOUBLE, root, tag, MPI_COMM_WORLD);
    MPI_Send(Mass, len[pId], MPI_DOUBLE, root, tag, MPI_COMM_WORLD);
  }
}

void Euler(double * Pos, double * Vel, const double * Mass, double * Acc, const double & dt, const int & N, function Accel, const int * len, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
  /*---------------------------------------------------------------------------
  Euler method to calculate position and velocity at next time step.
  -----------------------------------------------------------------------------
  Arguments:
    Pos   :   Position of bodies (1D vector).
    Vel   :   Velocity of bodies (1D vector).
    Mass  :   Mass of bodies (1D vector).
    Acc   :   Accelerations (1D vector).
    dt    :   Time step.
    N     :   Total number of bodies.
    Accel :   Function to calculate acceleration.
    len   :   Array with the number of bodies per node.
    tag   :   Message tag.
    pId   :   Process identity.
    nP    :   Number of processes.
    root  :   Root process.
    status:   Status object.
  ---------------------------------------------------------------------------*/
  Accel(Pos, Mass, Acc, len, N, tag, pId, nP, root, status);
  int ii, jj;
  for (ii = 0; ii < len[pId]; ii++){
    for (jj = 0; jj<3; jj++){
      Vel[3*ii+jj] += Acc[3*ii+jj]*dt;
      Pos[3*ii+jj] += Vel[3*ii+jj]*dt;
    }
  }
}

void PEFRL(double * Pos, double * Vel, const double * Mass, double *Acc, const double &dt, const int & N, function Accel, const int * len, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status){
  /*---------------------------------------------------------------------------
  Position Extended Forest-Ruth Like method to calculate position and velocity
  at next time step.
  -----------------------------------------------------------------------------
  Arguments:
    Pos   :   Position of bodies (1D vector).
    Vel   :   Velocity of bodies (1D vector).
    Mass  :   Mass of bodies (1D vector).
    Acc   :   Accelerations (1D vector).
    dt    :   Time step.
    N     :   Total number of bodies.
    Accel :   Function to calculate acceleration.
    len   :   Array with the number of bodies per node.
    tag   :   Message tag.
    pId   :   Process identity.
    nP    :   Number of processes.
    root  :   Root process.
    status:   Status object.
  ---------------------------------------------------------------------------*/
  int local_particles = len[pId];

  // Temporal arrays
  double * X = (double *) malloc(3*local_particles*sizeof(double));
  double * V = (double *) malloc(3*local_particles*sizeof(double));

  X = Pos;
  V = Vel;

  // PEFRL Parameters
  double xi = 0.1786178958448091e+0;
  double gamma = -0.2123418310626054e+0;
  double chi = -0.6626458266981849e-1;

  // Main loop
  int ii;
  for(ii = 0; ii < local_particles; ii++){
    X[3*ii] += xi*dt*V[3*ii];
    X[3*ii+1] += xi*dt*V[3*ii+1];
    X[3*ii+2] += xi*dt*V[3*ii+2];
  }
  
  Accel(X, Mass, Acc, len, N, tag, pId, nP, root, status);
  for (ii = 0; ii < local_particles; ii++){
    V[3*ii] +=  0.5*(1.-2*gamma)*dt*Acc[3*ii];
    V[3*ii+1] +=  0.5*(1.-2*gamma)*dt*Acc[3*ii+1];
    V[3*ii+2] +=  0.5*(1.-2*gamma)*dt*Acc[3*ii+2];
  }
  for (ii = 0; ii < local_particles; ii++){
    X[3*ii] += chi*dt*V[3*ii];
    X[3*ii+1] += chi*dt*V[3*ii+1];
    X[3*ii+2] += chi*dt*V[3*ii+2];
  }
  for (ii = 0; ii < 3*local_particles; ii++) Acc[ii] = 0.0;
  
  Accel(X, Mass, Acc, len, N, tag, pId, nP, root, status);
  for (ii = 0; ii < local_particles; ii++){
    V[3*ii] += gamma*dt*Acc[3*ii];
    V[3*ii+1] += gamma*dt*Acc[3*ii+1];
    V[3*ii+2] += gamma*dt*Acc[3*ii+2];
  }
  for (ii = 0; ii < local_particles; ii++){
    X[3*ii] += (1.-2*(chi+xi))*dt*V[3*ii];
    X[3*ii+1] += (1.-2*(chi+xi))*dt*V[3*ii+1];
    X[3*ii+2] += (1.-2*(chi+xi))*dt*V[3*ii+2];
  }
  for (ii = 0; ii < 3*local_particles; ii++) Acc[ii] = 0.0;

  Accel(X, Mass, Acc, len, N, tag, pId, nP, root, status);
  for (ii = 0; ii < local_particles; ii++){
    V[3*ii] += gamma*dt*Acc[3*ii];
    V[3*ii+1] += gamma*dt*Acc[3*ii+1];
    V[3*ii+2] += gamma*dt*Acc[3*ii+2];
  }
  for (ii = 0; ii < local_particles; ii++){
    X[3*ii] += chi*dt*V[3*ii];
    X[3*ii+1] += chi*dt*V[3*ii+1];
    X[3*ii+2] += chi*dt*V[3*ii+2];
  }
  for (ii = 0; ii < 3*local_particles; ii++) Acc[ii] = 0.0;
  
  Accel(X, Mass, Acc, len, N, tag, pId, nP, root, status);
  for (ii = 0; ii < local_particles; ii++){
    Vel[3*ii] = V[3*ii]+0.5*(1.-2*gamma)*dt*Acc[3*ii];
    Vel[3*ii+1] = V[3*ii+1]+0.5*(1.-2*gamma)*dt*Acc[3*ii+1];
    Vel[3*ii+2] = V[3*ii+2]+0.5*(1.-2*gamma)*dt*Acc[3*ii+2];
  }
  for (ii = 0; ii < local_particles; ii++){
    Pos[3*ii] = X[3*ii]+xi*dt*Vel[3*ii];
    Pos[3*ii+1] = X[3*ii+1]+xi*dt*Vel[3*ii+1];
    Pos[3*ii+2] = X[3*ii+2]+xi*dt*Vel[3*ii+2];
  }
}

void Evolution(std::ofstream & File, double *Pos, double * Vel, const double * Mass, double * Acc, const int * len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status, const int &steps, const double & dt, const int & jump, function Accel, Integrator evol){
  /*---------------------------------------------------------------------------
  Evolution of the system of bodies under gravitational interactions.
  -----------------------------------------------------------------------------
  Arguments:
    File  :   File where data is saved.
    Pos   :   Position of bodies (1D vector).
    Vel   :   Velocity of bodies (1D vector).
    Mass  :   Mass of bodies (1D vector).
    Acc   :   Accelerations (1D vector).
    len   :   Array with the number of bodies per node.
    N     :   Total number of bodies.
    tag   :   Message tag.
    pId   :   Process identity.
    nP    :   Number of processes.
    root  :   Root process.
    status:   Status object.
    steps :   Evolution steps.
    dt    :   Time step.
    jump  :   Data storage interval.
    Accel :   Function to calculate acceleration.
    evol  :   Integrator.
  ---------------------------------------------------------------------------*/
  
  if (pId==root){
    remove("Evolution.txt");
  }
  Save_data(File, Pos, Mass, len, N, tag, pId, nP, root, status);
  int ii, jj;
  for (ii = 0; ii < steps; ii++){
    evol(Pos, Vel, Mass, Acc, dt, N, Accel, len, tag, pId, nP, root, status);
    if ( ii%jump == 0) Save_data(File, Pos, Mass, len, N, tag, pId, nP, root, status);
    for (jj = 0; jj < 3*len[pId]; jj++) Acc[jj] = 0.0;
  }
}