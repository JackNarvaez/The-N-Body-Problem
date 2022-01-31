#include <iostream>
#include <cmath>
#include <math.h> 
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_double.h>


//g++ -std=c++17  fit.cpp -lgsl -lgslcblas

void LogNep(double*, int &n); 

int main(int argc, char **argv) {
  //Performs least-squares fit to a straight line
  double x[4] = {1,2,3,5};
  double y1[4] = {1,4,9,25};
  double y8[4] = {1,8,27,125};
  
  int n = sizeof(x)/sizeof(x[0]);
  
  LogNep(y1,n);
  LogNep(y8,n);
  LogNep(x,n);

  double c0, c1, cov00, cov01, cov11, sumsq;
  gsl_fit_linear (x,1,y1,1,n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

  double c08, c18, cov008, cov018, cov118, sumsq8;
  gsl_fit_linear (x,1,y8,1,n, &c08, &c18, &cov008, &cov018, &cov118, &sumsq8);
  
  double R1=gsl_stats_correlation (x,1,y1,1,n);
  double R8=gsl_stats_correlation (x,1,y8,1,n);
  
  std::cout<<"1 Process:\t"<<exp(c0)<<" x^"<<c1<<"\t"<<R1*R1<<"\n";
  std::cout<<"8 Processes:\t"<<exp(c08)<<" x^"<<c18<<"\t"<<R8*R8<<"\n";
  return 0;
}

void LogNep (double *vector, int &n){
  int i=0;
  for (i=0; i<n; i++){
    vector[i]=log(vector[i]);
  }
}
