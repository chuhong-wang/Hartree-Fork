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

        Matrix H_core(); 
        // S^(-1/2) = L_s * Λ_s^(-1/2) * L_s.T
        Matrix orthogonalization_matrix(const std::vector<double> &vec);
        Molecule get_molecule(); 

        void update_density(
            const Matrix& eigenVec_AObasis,
            int occupied_ao, 
            Matrix &orig_density); 

        // mat_out = orthogonalization_matrix_.transpose() * mat * orthogonalization_matrix_; 
        Matrix orthogonalize(const Matrix &mat, const Matrix &orthogonalization_matrix_); 
        std::pair<Matrix, Matrix> eigenVec_eigenVal(const Matrix &mat); 

        void update_Fork_matrix(const Matrix &h_core_, const Matrix &D, Matrix &orig_Fork); 

        double SCF_energy(const Matrix &density, const Matrix &h_core_, const Matrix &F_uv); 

    private:
        // convert a flattened-matrix vector back to a 2D Matrix
        Matrix vector_to_Matrix(const std::vector<double> &vec); 
       

        



        
};



