#include "include/scf.h" 
#include <sstream> 
#include <iostream> 
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

#define INDEX(i,j) (i>j) ? (i+1)*i/2+j : (j+1)*j/2+i 

// constructor of Scf initializes its molecule 
Scf::Scf(std::string geom_file, int q, int n_ao):mol(geom_file, q, n_ao) {
    mol.read_1e_integral("int1e_overlap"); 
    mol.read_1e_integral("int1e_kinetic"); 
    mol.read_1e_integral("int1e_nuc"); 
} 

Scf::~Scf() = default; 

// std::cout << "Core Hamiltonian" << std::endl; 
// for (auto i = 0; i!= scf.mol.nao_; ++i ){
//     for (auto j = 0; j!= scf.mol.nao_; ++j ){
//         auto ij = i*scf.mol.nao + j; 
//         std::printf("%8.5f \t", scf.T[ij] + scf.V[ij]); 
//     }
//     std::cout << std::endl; 
// } 



/**
 * diagonalization of overlap integral S 
 * SL_s = L_sΛ_s 
 * where L_s is the eigenvector as columns, 
 * and Λ_s is the eigenvalues on the diagonal of the matrix 
*/
Matrix Scf::orthogonalize_matrix(const Matrix &mat){
    int nao = mol.nao_; 
    Matrix mat_out(nao, nao);
    for (auto i = 0; i!=nao; ++i){
        for (auto j = 0; j!=nao; ++j){
            mat_out(i,j) = mol.one_e_overlap_en()[INDEX(i,j)]; 
        }
    }

    Eigen::SelfAdjointEigenSolver<Matrix> solver(mat_out); 
    Matrix evals = solver.eigenvalues(); 
    Matrix evecs = solver.eigenvectors(); 

    Matrix diag_evals(nao, nao); 
    for (auto i = 0;i!=nao; ++i) {
        diag_evals(i,i) = evals(i); 
    }
    Matrix matrix_sqrt = diag_evals.llt().matrixL(); // this is equivalent to squre root  
    diag_evals = matrix_sqrt.inverse(); 

    Matrix S_inv_sqrt = evecs * diag_evals * evecs.transpose(); 
    return S_inv_sqrt; 
}
