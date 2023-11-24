#pragma once 

#include <string> 
#include "constants.h"
#include "point.h"

class Molecule {
    private:
    int charge_;
    std::vector<int> zvals_; // atomic number of each atom
    int nao_; // number of atomic orbitals 
    std::vector<std::vector<double>> coords_;
    std::string point_group_ ; 

    public:
    // member functions 
    int num_atoms() {return coords_.size(); }
    void print_coord();
    double bond(int atom1, int atom2) const; 
    double angle(int atom1, int atom2, int atom3) const; 
    double torsion(int, int, int, int) const; 
    double out_of_plane_angle(int, int, int, int) const;
    double oop(int, int, int, int) const;  
    Point center_of_mass() const; 
    void translate(double, double, double); 
    void print_principal_mom_of_inertia() const; 

    // constructor
    Molecule() = default;  
    Molecule(std::string filename, int q, int n_ao=0); // read from file with user specified charge

    private:
    double m_unit_vec(int axis, int atom1, int atom2) const; 
    double* m_cross_prod(double* vec1, double* vec2) const; 
}; 


