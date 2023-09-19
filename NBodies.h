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

void read_data(const std::string &, double *, double *, double *);
void Gravitational_Acc(double *, const double *, const double *, const double *, const int &, const int &);
void Acceleration(const double *, const double *, double *, const int *, const int &, const int &, const int &, const int &, const int &, MPI_Status);
void Save_data(std::ofstream &, const double *, const double *, const int *, const int &, const int &, const int &, const int &, const int &, MPI_Status);
void Save_vec(std::ofstream &, const double *, const double *, const int &);
void PEFRL(double *, double *, const double *, double *, const double &, const int &, function, const int *, const int &, const int &, const int &, const int &, MPI_Status);
void Euler(double *, double *, const double *, double *, const double &, const int &, function, const int *, const int &, const int &, const int &, const int &, MPI_Status);
void Evolution(std::ofstream &, double *, double *, const double *, double *, const int *, const int &, const int &, const int &, const int &, const int &, MPI_Status, const int &, const double &, const int &, function, Integrator);

#endif // NBODIES_H_