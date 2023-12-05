#include "gtest/gtest.h"
#include "molecule.h"  
#include "scf.h"

TEST(MoleculeTest, BondTest) {
    // Create a Molecule instance for testing
    Molecule m("data/geom.dat", 0, "data/s.dat");

    // Perform tests on the bond function
    double expectedValue = 2.0786985874367461; 
    double tolerance = 1e-5; 

    EXPECT_NEAR(m.bond(0, 1), expectedValue, tolerance);
}

TEST(MoleculeTest, AngleTest) {
    // Create a Molecule instance for testing
    Molecule m("data/geom.dat", 0, "data/s.dat");

    // Perform tests on the angle function
    double expectedValue = 0.66322511575804732; 
    double tolerance = 1e-5; 
    EXPECT_NEAR(m.angle(0, 1, 2), expectedValue, tolerance);
    // Add more test cases as needed
}

TEST(MoleculeTest, ReadIntegralTest) {
    Molecule m("data/geom.dat", 0, "data/s.dat");
    m.read_1e_integral("int1e_overlap");
    double expectedValue = 1.0;
    double tolerance = 1e-5; 
    EXPECT_NEAR(m.one_e_overlap_en()[0], expectedValue, tolerance); 
    expectedValue = 0.236703936510848;
    EXPECT_NEAR(m.one_e_overlap_en()[1], expectedValue, tolerance); 
}

TEST(ScfTest, OrthogonizationTest){
    auto scf = Scf("data/geom.dat", 0, "data/s.dat"); 
    auto s_sqrt_int = scf.orthogonalization_matrix((scf.get_molecule()).one_e_overlap_en()); 
    double tolerance = 1e-5; 
    EXPECT_NEAR(s_sqrt_int(0,0), 1.0236346, tolerance);
    EXPECT_NEAR(s_sqrt_int(5,2), -0.1757583, tolerance);
}

TEST(ScfTest, HcoreTest){
    auto scf = Scf("data/geom.dat", 0, "data/s.dat"); 
    double tolerance = 1e-5; 
    auto h_core_ = scf.H_core(); 
    EXPECT_NEAR(h_core_(2,5), -1.6751501, tolerance); 
}

TEST(ScfTest, DensityTest){
    auto scf = Scf("data/geom.dat", 0, "data/s.dat"); 
    auto s_sqrt_int = scf.orthogonalization_matrix((scf.get_molecule()).one_e_overlap_en()); 
    
    double tolerance = 1e-5; 
    auto h_core_ = scf.H_core(); 
    auto F0_pr = scf.orthogonalize(h_core_, s_sqrt_int); 
    auto C0_pr = scf.eigenVec_eigenVal(F0_pr).first; 
    auto C0_ = s_sqrt_int*C0_pr; 

    auto mol = scf.get_molecule(); 
    Matrix D0(mol.get_nao(), mol.get_nao()); 
    scf.update_density(C0_, 5, D0); 
    EXPECT_NEAR(D0(0,0), 1.0650117, tolerance); 

    auto E0 = scf.SCF_energy(D0, h_core_, F0_pr); 
    EXPECT_NEAR(E0, -129.84539852990568, tolerance); 

    

}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
