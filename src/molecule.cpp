#include "molecule.h"
#include <fstream> 
#include <sstream>
#include <iostream>
#include <cmath> 
#include <cstdio>
#include <algorithm> 
#include <vector> 
#include <memory> 
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix; 

using namespace std; 

// constructor of Molecule from .dat file 
Molecule::Molecule(string filename, int q, int n_ao): charge_(q), nao_(n_ao){
    // read number of atoms and coordinates from file 
    ifstream fin(filename); 
    string line_zero;
    getline(fin, line_zero); 
    istringstream is(line_zero); 

    int natom; 
    is >> natom; 

    string line; 
    double zval, coord_0, coord_1, coord_2; 
    for(auto i=0; i!=natom;++i) {
        if (getline(fin, line)) {
            istringstream sin(line); 
            sin >>zval >> coord_0 >> coord_1 >> coord_2;   
            zvals_.push_back(zval); 
            coords_.push_back({coord_0, coord_1, coord_2}); 
        }
        else cerr << "inconsistent natom and coordinates\n"; 
    }
    fin.close(); 
}


// print coordinates of molecule (unit: Bohr)
void Molecule::print_coord(){
    cout << " coordinates of the molecule \n"; 
    for (auto i = 0; i!=coords_.size();++i) {
        printf ("%8.5f %8.5f %8.5f \n",coords_[i][0], coords_[i][1], coords_[i][2]); 
    }
}

// Params:
//  arg1: index of the first atom
//  arg2: index of the second atom
// Returns: bond length = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
double Molecule::bond(int atom1, int atom2) const{
    if (atom1>=coords_.size() || atom2>=coords_.size()){
        cerr<< "invalid atom index" << endl; 
    } 
    double sum_sqr = 0;
    for (auto i=0; i!=3; ++i) {
        auto diff_i = coords_[atom1][i] - coords_[atom2][i];
        sum_sqr += diff_i * diff_i; 
    }
    return sqrt(sum_sqr); 
}

double Molecule::m_unit_vec(int axis, int atom1, int atom2) const{
    return (coords_[atom2][axis] - coords_[atom1][axis])/bond(atom1, atom2); 
}

// Params: indices of three atoms i-j-k
// Returns: angle formed by i-j-k
double Molecule::angle(int atom1, int atom2, int atom3) const{
    return acos(m_unit_vec(0, atom2, atom1) * m_unit_vec(0, atom2, atom3) +
                m_unit_vec(1, atom2, atom1) * m_unit_vec(1, atom2, atom3) +
                m_unit_vec(2, atom2, atom1) * m_unit_vec(2, atom2, atom3)); 
}

double* Molecule::m_cross_prod(double* vec1, double* vec2) const{
    double* ret = new double[3]; 
    ret[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1]; 
    ret[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2]; 
    ret[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
    return ret; 
}

// theta_ijkl = arcsin( (unit_vec(kj) X unit_vec(kl)/ sin(angle(j,k,l)) dot unit_vec(k, i)))
double Molecule::out_of_plane_angle(int i, int j, int k, int l) const {
    double vec_kj[3] = {m_unit_vec(0, k, j), m_unit_vec(1, k, j), m_unit_vec(2, k, j)}; 
    double vec_kl[3] = {m_unit_vec(0, k, l), m_unit_vec(1, k, l), m_unit_vec(2, k, l)};
    auto cross_prod = m_cross_prod(vec_kj, vec_kl); 
    auto numerator = cross_prod[0]*m_unit_vec(0, k, i) + 
                    cross_prod[1]*m_unit_vec(1, k, i) + 
                    cross_prod[2]*m_unit_vec(2, k, i); 
    delete cross_prod; 
    auto theta = numerator / sin(angle(j,k,l)); 

    // handles the value outside [-1, 1] caused by numerical precision 
    if (theta > 1.0) {theta = asin(1.0); }
    else if (theta < -1.0) {theta = asin(-1.0);}
    else {theta = asin(theta); }
    return theta; 
}

double Molecule::torsion(int i, int j, int k, int l) const{
    // compute the angle b/w a-b-c and b-c-d
    double vec_ij[3] = {m_unit_vec(0, i, j), m_unit_vec(1, i, j), m_unit_vec(2, i, j)}; 
    double vec_jk[3] = {m_unit_vec(0, j, k), m_unit_vec(1, j, k), m_unit_vec(2, j, k)}; 
    double vec_kl[3] = {m_unit_vec(0, k, l), m_unit_vec(1, k, l), m_unit_vec(2, k, l)}; 
    auto cross_prod_ij_jk = m_cross_prod(vec_ij, vec_jk); 
    auto cross_prod_jk_kl = m_cross_prod(vec_jk, vec_kl);
    auto dot_prod = cross_prod_ij_jk[0] * cross_prod_jk_kl[0] + 
                    cross_prod_ij_jk[1] * cross_prod_jk_kl[1] + 
                    cross_prod_ij_jk[2] * cross_prod_jk_kl[2]; 


    auto denom = sin(angle(i,j,k))* sin(angle(j,k,l)); 
    auto tau = dot_prod/denom; 
    if (tau > 1.0) {tau = acos(1.0); }
    else if (tau < -1.0) {tau = acos(-1.0); }
    else {tau = acos(tau); }
    
    // compute the sign of torsion: clockwise->positive, counterclock->negative 
    // sin(tau) >0 -> positive, <0 -> negative 
    int sign = 1; 
    auto cross_prod_3 = m_cross_prod(cross_prod_ij_jk, cross_prod_jk_kl); 
    auto dot_prod_sin = vec_jk[0] * cross_prod_3[0] + 
                        vec_jk[1] * cross_prod_3[1] +
                        vec_jk[2] * cross_prod_3[2]; 
    if (dot_prod_sin < 0) sign = -1; 
    delete cross_prod_3; 
    delete cross_prod_ij_jk;
    delete cross_prod_jk_kl; 

    return tau*sign; 
}

Point Molecule::center_of_mass() const{
    vector<double> com(3, 0.0); 
    double sum_atomic_masses = 0.0; 
    for (auto i = 0; i!= coords_.size(); ++i){
        com[0] += atomic_masses[zvals_[i]] * coords_[i][0]; 
        com[1] += atomic_masses[zvals_[i]] * coords_[i][1]; 
        com[2] += atomic_masses[zvals_[i]] * coords_[i][2]; 
        sum_atomic_masses += atomic_masses[zvals_[i]]; 
    }
    com[0]/=sum_atomic_masses; 
    com[1]/=sum_atomic_masses; 
    com[2]/=sum_atomic_masses; 

    return Point(com); 
}

void Molecule::translate(double a, double b, double c){
    for (auto i =0; i!=coords_.size(); ++i){
        coords_[i][0]+= a;
        coords_[i][1]+= b;
        coords_[i][2]+= c; 
    }
}

void Molecule::print_principal_mom_of_inertia() const{
    // calculate mom of inertia
    Matrix I(3,3); 
    for (auto i = 0; i!=coords_.size(); ++i) {
        I(0,0) += atomic_masses[zvals_[i]]*(coords_[i][1]* coords_[i][1] + coords_[i][2]* coords_[i][2]); 
        I(1,1) += atomic_masses[zvals_[i]]*(coords_[i][0]* coords_[i][0] + coords_[i][2]* coords_[i][2]);
        I(2,2) += atomic_masses[zvals_[i]]*(coords_[i][0]* coords_[i][0] + coords_[i][1]* coords_[i][1]);
        I(0,1) += atomic_masses[zvals_[i]]*coords_[i][0]* coords_[i][1]; 
        I(1,2) += atomic_masses[zvals_[i]]*coords_[i][1]* coords_[i][2]; 
        I(0,2) += atomic_masses[zvals_[i]]*coords_[i][0]* coords_[i][2];
    }
    I(0,1) *= -1; 
    I(1,0) = I(0,1);
    I(0,2) *= -1;
    I(2,0) = I(0,2); 
    I(1,2) *= -1; 
    I(2,1) = I(1,2); 
    cout << "moment of inertia: \n"; 
    cout << I << endl; 
    Eigen::SelfAdjointEigenSolver<Matrix> solver(I); 
    Matrix evecs = solver.eigenvectors(); 
    Matrix evals = solver.eigenvalues(); 
    cout << "principal moment of inertia (amu * bohr^2): \n" << evals << endl; 
}