#include "include/scf.h" 
#include <sstream> 
#include <iostream> 
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

// constructor of Scf initializes its molecule 
Scf::Scf(string geom_file, int q, int n_ao):mol(geom_file, q, n_ao) {} 

Scf::~Scf() = default; 

void Scf::read_energy_scalar(string filename){
    std::ifstream fin(filename); 
    string s; 
    std::getline(fin, s); 
    std::istringstream iss(s); 
    iss >> neuc;
}


void Scf::read_1e_integral(string filename, vector<double>& vec){
    std::ifstream fin(filename); 
    string line; 
    int idx1, idx2; 
    double val; 
    while (getline(fin, line)){
        std::istringstream iss(line); 
        iss >> idx1 >> idx2 >> val; 
        int idx_flt = idx1*mol.nao + idx2; 
        vec[idx_flt] = val; 
    }
}



/**
 * 2-D symmetrical matrix can be saved as a 1-D array with reduced space 
 *      (0,0)                           0
 *      (1,0) (1,1)                     1   2
 *      (2,0) (2,1) (2,2)          =>   3   4   5
 *      (3,0) (3,1) (3,2) (3,3)         6   7   8   9 
 *                 ...                          ...
 * In this way, a N*N matrix is reduced to (N+1)*N/2 units of memory 
 * 
 * Similarly, ijkl with 8-order symmetry can be reduced by treating ij as dim_0 and kl and dim_1
*/
void Scf::twoElectronIntegral(string filename){
    int max_ij = mol.nao*(mol.nao-1)/2 + (mol.nao-1); 
    int eri_count = (max_ij+1)*max_ij/2 + max_ij; 
    ee = std::vector<double>(eri_count+1, 0);  // e.g. (0,0,0,0) -> (7,7,7,7) extra 1 slot for index 0

    std::ifstream fin1(filename); 
    string line; 
    int i,j,k,l; 
    double val; 
    while (getline(fin1, line)){
        std::istringstream iss(line); 
        iss >> i >> j >> k >> l >> val; 
        --i; --j; --k; --l; 
        int ij = i>j? (i+1)*i/2 + j : (j+1)*j/2+i; 
        int kl = k>l? (k+1)*k/2 + l : (l+1)*l/2+k; 
        int ijkl = ij>kl? (ij+1)*ij/2 + kl : (kl+1)*kl/2 + ij; 
        // std::cout << ij << " " << kl << " " << ijkl << std::endl; 
        ee[ijkl] = val; 
        // std::cout << val << std::endl; 
    }
}