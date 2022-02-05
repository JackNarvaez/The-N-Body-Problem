#ifndef NBODIES_H_
#define NBODIES_H_
#include <cassert>

void read_NParticles(const std::string &File_address, int& N);
void read_data(const std::string &File_address, std::vector<double>& Pos, std::vector<double>& Mom);
void Initial_state(const std::string &File_address, std::vector<double>& Pos, std::vector<double>&Mom, std::vector<int>&len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Gravitational_force(std::vector<double>& Force, const std::vector<double>& vec0, const std::vector<double>& vec1, const int & len0, const int & len1);
void Total_Force(const std::vector<double>& Pos, std::vector<double>& Force, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Save_data(std::ofstream& File, const std::vector<double>& Pos, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Save_vec(std::ofstream& File, const std::vector<double>& Vec, const int & len);
void Euler(std::vector<double>& Pos, std::vector<double>& Mom, const std::vector<double>& Force, const double &dt, const int & N);
void Evolution(std::ofstream& File, std::vector<double>& Pos, std::vector<double>& Mom, std::vector<double>& Force, const std::vector<int>& len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status, const int &steps, const double & dt, const int & jump);

#endif // NBODIES_H_