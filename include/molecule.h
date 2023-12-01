#pragma once 

#include <string> 
#include <vector> 
#include "constants.h"
#include "point.h"

class Molecule {
    friend class Scf; 
    private:
    int charge_;
    std::vector<int> zvals_; // atomic charge of each atom
    int nao_; // number of atomic orbitals 
    std::vector<std::vector<double>> coords_;
    std::string point_group_ ; 

    double enuc_;  
    std::vector<double> S_uv_; 
    std::vector<double> T_uv_;
    std::vector<double> V_uv_; 

    std::vector<double> eri_; 

    public:
    // member functions 
    int num_atoms() {return coords_.size(); }
    int get_nao() {return nao_; }
    void print_coord();
    double bond(int atom1, int atom2) const; 
    double angle(int atom1, int atom2, int atom3) const; 
    double torsion(int, int, int, int) const; 
    double out_of_plane_angle(int, int, int, int) const;
    double oop(int, int, int, int) const;  
    Point center_of_mass() const; 
    void translate(double, double, double); 
    void print_principal_mom_of_inertia() const; 

    double nuclear_repulsion_en() const; 
    std::vector<double> one_e_overlap_en() const; 
    std::vector<double> one_e_kinetic_en() const; 
    std::vector<double> one_e_nuclearAttraction_en() const; 

    std::vector<double> two_e_integral() const; 

    // constructor
    Molecule() = default;  
    Molecule(std::string filename, int q, int n_ao=0); // read from file with user specified charge
    void read_1e_integral(std::string integral_type); 
    void read_2e_integral(std::string filename); 
    void read_energy_scalar(std::string filename);  // read nuclear repulsion energy 

    private:
    double m_unit_vec(int axis, int atom1, int atom2) const; 
    double* m_cross_prod(double* vec1, double* vec2) const; 
    void file_to_vector(std::string filename, std::vector<double>& vec); 
    

}; 


