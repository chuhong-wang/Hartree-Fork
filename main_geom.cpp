#include "include/molecule.h"
#include <cstdio> 
int main() {
    Molecule m("data/geom.dat", 0); 
    m.print_coord(); 
    printf("interatomic bonds \n"); 
    for (auto i=0; i!=m.num_atoms(); ++i){
        for (auto j = 0; j!=i; ++j) {
            printf("%d %d %8.5f\n", i,j, m.bond(i,j));
        }
    }

    printf("bond angle \n"); 
    for (auto i=0; i!=m.num_atoms();++i) {
        for (auto j=0; j!=i; ++j) {
            for (auto k=0; k!=j; ++k){
                if(m.bond(i,j) < 4.0 && m.bond(j,k) < 4.0) {
                    printf("%d %d %d %8.5f \n", i, j, k, m.angle(i, j, k)/acos(-1)*180); 
                }
                
            }
        }
    }

    printf("out-of-plane angle \n"); 
    for (auto i=0; i!=m.num_atoms();++i) {
        for (auto j=0; j!=m.num_atoms(); ++j) {
            for (auto k=0; k!=m.num_atoms(); ++k){
                for (auto l=0; l!=j; ++l) {
                    if (i!=j && i!=k && i!=l && j!=k && k!=l  
                        && m.bond(i,k) < 4.0 && m.bond(k,j) < 4.0 && m.bond(k,l) < 4.0){
                        printf("%d %d %d %d %8.5f \n", 
                        i, j, k, l, m.out_of_plane_angle(i, j, k, l)/acos(-1)*180);
                    }
                }    
            }
        }
    }

    printf("Torsional angles\n");
    for(int i=0; i < m.num_atoms(); i++) {
        for(int j=0; j < i; j++) {
        for(int k=0; k < j; k++) {
            for(int l=0; l < k; l++) {
            if(m.bond(i,j) < 4.0 && m.bond(j,k) < 4.0 && m.bond(k,l) < 4.0)
                printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, m.torsion(i,j,k,l)*(180.0/acos(-1.0)));
            }
        }
        }
    }

    printf("Center of Mass\n");
    auto com = m.center_of_mass(); 
    printf("%8.5f %8.5f %8.5f\n", com.get_x(), com.get_y(), com.get_z());
    
    m.translate(-com.get_x(), -com.get_y(), -com.get_z()); 
    m.print_principal_mom_of_inertia(); 
    return 0; 
}