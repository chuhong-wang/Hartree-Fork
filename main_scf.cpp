#include "include/scf.h"

#include <vector>
#include <string>
#include <iostream>
#include <set> 
#include <cstdio> 

#define INDEX(i,j) (i>j) ? (i+1)*i/2+j : (j+1)*j/2+i 

/**
 * ij = (i+1)*i/2 + j
 * (i+1)*i gets calculated multiple times -> 
 * we can also store its value in an pre-computed array for fast lookup
 * 
*/

int main() {

    auto scf = Scf("data/geom.dat", 0, 7); 

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
     * SL_s = L_sΛ_s where L_s is the eigenvector and Λ_s is the eigenvalue
    */
    int nao = scf.mol.nao; 
    Matrix I(nao, nao); //re-name 
    for (auto i = 0; i!=nao; ++i){
        for (auto j = 0; j!=nao; ++j){
            I(i,j) = scf.overlap[i][j];
        }
    }
    Eigen::SelfAdjointEigenSolver<Matrix> solver(I); 
    Matrix evals = solver.eigenvalues(); 
    std::cout << evals << std::endl; 
    Matrix evecs = solver.eigenvectors(); 

    Matrix diag_evals(nao, nao); 
    for (auto i = 0;i!=nao; ++i) {
        diag_evals(i,i) = evals(i); 
    }
    Matrix matrix_sqrt = diag_evals.llt().matrixL(); // this is equivalent to squre root  
    diag_evals = matrix_sqrt.inverse(); 
    std::cout << "diagonal matrix of eigenvalues" << std::endl; 
    std::cout << diag_evals << std::endl; 
    Matrix S_inv_sqrt = evecs * diag_evals * evecs.transpose(); 
    std::cout << "\n The S^(-1/2) matrix :"<<std::endl; 
    std::cout << S_inv_sqrt << std::endl; 

    /**
     * 1. initial (guess) Fork matrix 
     *    F_0' = S^(-1/2).T * H_core * S^(-1/2)
     * 2. diagonalize the Fork matrix
     *    F_0' C_0' = C_0' ε_0'
     * 3. transform the eigenvectors into the original (non-orthogonal) AO basis
     *    C_0 = S^(-1/2)C_0'
     * 4. build the density matrix using the occupied MOs 
     * 
    */
    Matrix H_core(nao, nao); 
    for (auto i = 0; i!=nao; ++i){
        for (auto j = 0; j!=nao; ++j){
            H_core(i,j) = scf.kineticE[i][j] + scf.attraction[i][j];
        }
    }
    std::cout << "H_core" << std::endl; 
    std::cout << H_core << std::endl; 
    Matrix F0_pr = S_inv_sqrt.transpose() * H_core * S_inv_sqrt; 
    std::cout << "Transformed Fork Matrix (orthogonal basis)"<<std::endl; 
    std::cout << F0_pr << std::endl; 
    solver.compute(F0_pr); 
    evecs = solver.eigenvectors();
    evals = solver.eigenvalues(); 
    Matrix C0 = S_inv_sqrt * evecs; 
    std::cout << "initial MO coefficient (AO basis)"<< std::endl; 
    std::cout << C0 << std::endl; // TODO: some signs don't match with the answer

    Matrix C0_occAO = C0.block(0, 0, 7, 5); 
    Matrix D = C0_occAO * C0_occAO.transpose(); 
    std::cout << "initial guess of density matrix" << std::endl; 
    std::cout << D << std::endl; 

    // initial guess: D=0
    for (auto i =0; i!=nao; ++i){
        for (auto j = 0; j!=nao; ++j){
            D(i,j)=0; 
        }
    } 

    // initial SCF energy 
    auto E0 = (D.array()*(H_core + F0_pr).array()).sum(); 
    std::cout <<"E0_elec: " << E0 <<" E_nuc: "<< scf.neuc <<std::endl; 

    // compute new Fork Matrix 
    std::vector<std::vector<double>> F(nao, std::vector<double>(nao, 0)); 
    Matrix F0_mat(nao, nao); 
    for (int i = 0; i!=nao; ++i){
        for (int j = 0; j!=nao; ++j){
            F[i][j] = H_core(i, j); 
            for (int k=0; k!=nao; ++k){
                for (int l=0; l!=nao; ++l){
                    int ij = INDEX(i,j); 
                    int kl = INDEX(k,l); 
                    int ijkl = INDEX(ij, kl);                        
                    int ik = INDEX(i, k); 
                    int jl = INDEX(j, l); 
                    int ikjl = INDEX(ik, jl); 
                    F[i][j] += D(k,l) * (2.0*scf.ee[ijkl] - scf.ee[ikjl]); 
                    F0_mat(i,j) = F[i][j]; 
                }
            }
        }
    }

    std::cout << "New Fork Matrix: "<<std::endl; 
    for (auto i =0; i!=nao ;++i){
        for (auto j=0; j!=nao; ++j) {
            std::printf("%8.5f ", F[i][j]); 
        }
        std::cout << std::endl; 
    }
    for (auto iter = 0; iter !=10; iter++){
        F0_pr = S_inv_sqrt.transpose() * F0_mat * S_inv_sqrt; 
        solver.compute(F0_pr); 
        evecs = solver.eigenvectors();
        evals = solver.eigenvalues(); 
        C0 = S_inv_sqrt * evecs; 
        C0_occAO = C0.block(0, 0, 7, 5); 
        D = C0_occAO * C0_occAO.transpose(); 
        E0 = (D.array()*(H_core + F0_pr).array()).sum(); 
        std::cout <<"E0_elec: " << E0 <<" E_nuc: "<< scf.neuc <<std::endl; 
        
        for (int i = 0; i!=nao; ++i){
            for (int j = 0; j!=nao; ++j){
                F[i][j] = H_core(i, j); 
                for (int k=0; k!=nao; ++k){
                    for (int l=0; l!=nao; ++l){
                        int ij = INDEX(i,j); 
                        int kl = INDEX(k,l); 
                        int ijkl = INDEX(ij, kl);                        
                        int ik = INDEX(i, k); 
                        int jl = INDEX(j, l); 
                        int ikjl = INDEX(ik, jl); 
                        F[i][j] += D(k,l) * (2.0*scf.ee[ijkl] - scf.ee[ikjl]); 
                        F0_mat(i,j) = F[i][j]; 
                    }
                }
            }
        }
        // std::cout << "New Fork Matrix: "<<std::endl; 
        // for (auto i =0; i!=nao ;++i){
        //     for (auto j=0; j!=nao; ++j) {
        //         std::printf("%8.5f ", F[i][j]); 
        //     }
        //     std::cout << std::endl; 
        // }

    }
    
    return 0; 

}