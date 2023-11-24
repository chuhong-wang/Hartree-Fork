#ifndef SCF_H
#define SCF_H

#include "molecule.h"
#include <fstream> 
#include <string> 
#include <vector> 
#include <memory> 

using std::string, std::vector; 

class Scf {
    public:
        Molecule mol; 
        double neuc; 
        vector<double> S; 
        vector<double> T;
        vector<double> V; 
        vector<double> eri; 


        // constructor 
        Scf() = default; 
        Scf(string geom_file, int q, int n_ao); // charge and number of AO of the molecule
        
        // desctructor 
        ~Scf(); 

        // functions for reading input data 
        void read_energy_scalar(string filename = "data/enuc.dat"); // TODO: remove default val
        // one-electron integrals 
        void read_1e_integral(string filename, vector<double>& out_data); // 
        // two-electron integrals 
        void twoElectronIntegral(string filename, vector<double>&); 

};



#endif 