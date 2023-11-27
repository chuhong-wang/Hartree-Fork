#include "scf.h"

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

    auto s_sqrt_int = scf.orthogonalization_matrix((scf.get_molecule()).one_e_overlap_en()); 
    auto h_core_ = scf.H_core(); 
    auto F0_pr = scf.orthogonalize(h_core_, s_sqrt_int); 
    auto C0_pr = scf.eigenVec_eigenVal(F0_pr).first; 
    auto C0_ = s_sqrt_int*C0_pr; 

    Matrix D0(7,7); 
    scf.update_density(C0_, 5, D0); 
    double E0 = scf.SCF_energy(D0, h_core_, F0_pr); 
    std::cout << E0 << std::endl; 

    for(auto i = 0; i!=20; ++i){
        scf.update_Fork_matrix(h_core_, D0, F0_pr); 
        scf.update_density(s_sqrt_int*scf.eigenVec_eigenVal(F0_pr).first, 5, D0); 
        double E = scf.SCF_energy(D0, h_core_, F0_pr);  
        std::cout << E << std::endl;
    }
    return 0; 

}