#include "scf.h"

#include <vector>
#include <string>
#include <iostream>
#include <set> 
#include <cstdio> 
#include <algorithm> 

#define INDEX(i,j) (i>j) ? (i+1)*i/2+j : (j+1)*j/2+i 

int main() {

    auto scf = Scf("data/geom.dat", 0, "data/s.dat"); 

    int nao = scf.get_molecule().get_nao(); 
    int ndocc = 5; 

    auto s_sqrt_int = scf.orthogonalization_matrix((scf.get_molecule()).one_e_overlap_en()); 
    auto h_core_ = scf.H_core(); 
    auto F0_pr = scf.orthogonalize(h_core_, s_sqrt_int); 
    auto C0_pr = scf.eigenVec_eigenVal(F0_pr).first; 
    auto C0_ = s_sqrt_int*C0_pr; 

    Matrix D0(nao,nao); 
    scf.update_density(C0_, ndocc, D0); 
    double E0 = scf.SCF_energy(D0, h_core_, F0_pr); 
    std::cout << E0 << std::endl; 
    std::cout << "initial density matrix " << std::endl; 
    std::cout << D0 << std::endl; 
    for(auto i = 0; i!=20; ++i){
        scf.update_Fork_matrix(h_core_, D0, F0_pr); 
        scf.update_density(s_sqrt_int*scf.eigenVec_eigenVal(F0_pr).first, 5, D0); 
        // std::cout << D0 << std::endl; 
        double E = scf.SCF_energy(D0, h_core_, F0_pr);  
        std::cout << E << std::endl;
    }
    std::cout << "density matrix " << std::endl; 
    std::cout << D0 << std::endl; 
    auto eps = scf.eigenVec_eigenVal(F0_pr).second; 
    auto eri_MO_ = scf.Noddy_algo(C0_); 
    double Emp2 = 0.0;
    int i, j, a, b, ia, jb, ib, ja, iajb, ibja; 
    for(i=0; i < ndocc; i++) {
        for(a=ndocc; a < nao; a++) {
        ia = INDEX(i,a);
        for(j=0; j < ndocc; j++) {
            ja = INDEX(j,a);
            for(b=ndocc; b < nao; b++) {
            jb = INDEX(j,b);
            ib = INDEX(i,b);
            iajb = INDEX(ia,jb);
            ibja = INDEX(ib,ja);
            Emp2 += scf.get_molecule().two_e_integral()[iajb] * 
                    (2 * scf.get_molecule().two_e_integral()[iajb] - scf.get_molecule().two_e_integral()[ibja])/
                    (eps(i) + eps(j) - eps(a) - eps(b));
            }
        }
        }
    }
    std::cout << Emp2 << std::endl; 

    return 0; 

}