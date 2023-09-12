#ifndef NBODIES_H_
#define NBODIES_H_
#include <cassert>

struct _body {
    double *m;
    double *r;
    double *v;
    double *a;
};

typedef struct _body body;

using function = void(const double *, const double *, double *, const int *, const int &, const int &, const int &, const int &, const int &, MPI_Status);
using Integrator = void(double *, double *, const double *, double *, const double &, const int &, function, const int *, const int &, const int &, const int &, const int &, MPI_Status);

void read_data(const std::string &File_address, double * Pos, double * Vel, double * Mass);
void Gravitational_Acc(double * Acc, const double * Pos0, const double * Pos1, const double * Mass1, const int & len0, const int & len1);
void Acceleration(const double * Pos, const double * Mass, double * Acc, const int* len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Save_data(std::ofstream& File, const double * Pos, const int * len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Save_vec(std::ofstream& File, const double * Vec, const int & len);
void PEFRL(double * Pos, double * Vel, const double * Mass, double * Acc, const double &dt, const int & N, function Accel, const int * len, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Euler(double * Pos, double * Vel, const double * Mass, double * Acc, const double &dt, const int & N, function Accel, const int * len, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status);
void Evolution(std::ofstream& File, double * Pos, double * Vel, const double * Mass, double * Acc, const int * len, const int & N, const int & tag, const int & pId, const int & nP, const int & root, MPI_Status status, const int &steps, const double & dt, const int & jump, function Accel, Integrator evol);

#endif // NBODIES_H_