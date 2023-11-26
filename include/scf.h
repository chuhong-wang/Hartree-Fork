#pragma once 

#include "molecule.h"
#include <fstream> 
#include <string> 
#include <vector> 
#include <memory> 

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix; 

class Scf {
    private:
        Molecule mol; 
    public:
        // constructor 
        Scf() = default; 
        Scf(std::string geom_file, int q, int n_ao); // charge and number of AO of the molecule

        Matrix orthogonalize_matrix(const Matrix &mat);
        
    private:
        Matrix vector_to_Matrix(const std::vector<double> &vec); 
        


};



