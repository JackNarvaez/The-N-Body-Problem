#ifndef DIFFUSION_H_
#define DIFFUSION_H_
#include <cassert>

void read_data(const std::string &File_address, std::vector<double>& Pos, std::vector<double>& Mom);
void Initial_state(const std::string &File_address, std::vector<double>& Pos, std::vector<double>&Mom, std::vector<int>&len, int N, int tag, int pId, int nP, int root, MPI_Status status);
void Total_Force(std::vector<double>& Pos, std::vector<double>& Force, std::vector<int>& len, int N, int tag, int pId, int nP, int root, MPI_Status status);
void Gravitational_force(std::vector<double>& Force, const std::vector<double>& vec0, const std::vector<double>& vec1, int len0, int len1);
void Save_data(std::ofstream& File, std::vector<double>& Force, const std::vector<int>& len, int N, int tag, int pId, int nP, int root, MPI_Status status);
void Save_vec(std::ofstream& File, std::vector<double>& Vec, int len);

#endif // DIFFUSION_H_