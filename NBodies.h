#ifndef NBODIES_H_
#define NBODIES_H_
#include <cassert>

using function = void(const std::vector<double>&, const std::vector<double>&, std::vector<double>&, const std::vector<int>&, const int &, const int &, const int &, const int &, const int &, MPI_Status);
using Integrator = void(std::vector<double>&, std::vector<double>&, const std::vector<double>&, std::vector<double>&, const double &, const int &, function, const std::vector<int>&, const int &, const int &, const int &, const int &, MPI_Status);

void read_data(const std::string &File_address, std::vector<double>& Pos, std::vector<double>& Vel, std::vector<double>& Mass);
void Gravitational_Acc(std::vector<double>& Acc, const std::vector<double>& Pos0, const std::vector<double>& Pos1, const std::vector<double>& Mass1, const int & len0, const int & len1);
void Acceleration(const std::vector<double>& Pos, const std::vector<double>& Mass, std::vector<double>& Acc, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Save_data(std::ofstream& File, const std::vector<double>& Pos, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Save_vec(std::ofstream& File, const std::vector<double>& Vec, const int & len);
void Euler(std::vector<double>& Pos, std::vector<double>& Vel, const std::vector<double>& Mass, std::vector<double>& Acc, const double &dt, const int & N, function Accel, const std::vector<int>& len, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void PEFRL(std::vector<double>& Pos, std::vector<double>& Vel, const std::vector<double>& Mass, std::vector<double>& Acc, const double &dt, const int & N, function Accel, const std::vector<int>& len, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Evolution(std::ofstream& File, std::vector<double>& Pos, std::vector<double>& Vel, const std::vector<double>& Mass, std::vector<double>& Acc, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status, const int &steps, const double & dt, const int & jump, function Accel, Integrator evol);

#endif // NBODIES_H_