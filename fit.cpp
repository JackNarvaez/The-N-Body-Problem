#include <iostream>
#include <cmath>
#include <fstream> 
#include <sstream> 
#include <vector>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_double.h>


//g++ -std=c++17  fit.cpp -lgsl -lgslcblas

void read_data(const std::string &File_address, std::vector<double>& vector, int column);
void Map(const std::vector<double>& x, std::vector<double>& y, const int &n);

int main(int argc, char **argv){
  std::vector<double> x;
  std::vector<double> y;
  
  std::string ext{"strong"}, format{".txt"}; 
  std::string data = ext+argv[1]+format;
  read_data(data, x, 3);
  read_data(data, y, 2);
  int n = x.size();

  Map(x, y, n);

  double c0, c1, cov00, cov01, cov11, sumsq;
  gsl_fit_linear (x.data(),1,y.data(),1,n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
  
  double R1 = gsl_stats_correlation (x.data(),1,y.data(),1,n);
  
  std::cout<<c0 << "x + " << c1 <<"x^2"<<"\t"<<R1*R1<<"\n";
  return 0;
}

void read_data(const std::string &File_address, std::vector<double>& vector, int column){
  /*---------------------------------------------------------------------------
  read_data:
  Read the column <column> from <File_address>.
  -----------------------------------------------------------------------------
  Arguments:
    File_address:   File address from which the data is read.
    Vector      :   Vector to store the data.
    column      :   Number of column for reading.
  ---------------------------------------------------------------------------*/
  std::ifstream File;
  File.open (File_address, std::ifstream::in);    // Open file
  std::string line;
  std::istringstream data;
	while (!File.eof()){
	std::getline(File,line);
  // Omit empty lines and comments
	if (line.length() == 0 || line[0] == '#'){
		continue;
  }else{
    std::istringstream iss(line);   // Separate line in columns
    std::string data;        
    for (int ii=0; ii < column; ii++){
      iss>>data;
    }
    vector.push_back(atof(data.c_str()));
  }
  }
  File.close();
}

void Map(const std::vector<double>& x, std::vector<double>& y, const int &n){
  /*---------------------------------------------------------------------------
  Map:
  Map y(x) to z(x)=y(x)/x .
  -----------------------------------------------------------------------------
  Arguments:
    x:   Domain.
    y:   y(x).
    n:   Size
  ---------------------------------------------------------------------------*/
  for (int ii=0; ii<n; ii++){
    y[ii]/=x[ii];
  }
}