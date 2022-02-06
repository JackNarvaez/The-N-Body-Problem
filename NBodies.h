#ifndef NBODIES_H_
#define NBODIES_H_
#include <cassert>

void read_NParticles(const std::string &File_address, int& N);
void read_data(const std::string &File_address, std::vector<double>& Pos, std::vector<double>& Vel);
void Initial_state(const std::string &File_address, std::vector<double>& Pos, std::vector<double>& Vel, std::vector<int>&len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Gravitational_Acc(std::vector<double>& Acc, const std::vector<double>& vec0, const std::vector<double>& vec1, const int & len0, const int & len1);
void Acceleration(const std::vector<double>& Pos, std::vector<double>& Acc, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Save_data(std::ofstream& File, const std::vector<double>& Pos, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Save_vec(std::ofstream& File, const std::vector<double>& Vec, const int & len);
void Euler(std::vector<double>& Pos, std::vector<double>& Vel, const std::vector<double>& Acc, const double &dt, const int & N);
void Evolution(std::ofstream& File, std::vector<double>& Pos, std::vector<double>& Vel, std::vector<double>& Acc, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status, const int &steps, const double & dt, const int & jump);

#endif // NBODIES_H_