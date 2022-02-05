/*---------------------------------------------
Create a file with the initial configuration of
a system of particles enclosed in a box of unit
side, the mass is constrained between 0 and 1.
---------------------------------------------*/

#include<iostream>
#include <fstream>  // std::ofstream
#include <random>

int main(int argc, char **argv){
    std::ofstream File;
    double seed{0};
    int particles = std::atoi(argv[1]);  // Number of particles
    File.open ("Random.txt", std::ofstream::out); // Store initial state
    
    File << "#\t" << particles << std::endl;    // Write in first line the number of particles.

    std::mt19937 g1(seed);
    std::uniform_real_distribution<double> distr(0.0,1.0);

    for(int ii=0; ii < particles; ii++){

        for (int jj=0; jj < 7; jj++){
            File << distr(g1) << "\t";
        }
        File << std::endl;
    }
    File.close();

    return 0;
}
