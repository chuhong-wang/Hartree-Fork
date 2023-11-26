#include "gtest/gtest.h"
#include "molecule.h"  // Include the header file for your Molecule class

TEST(MoleculeTest, BondTest) {
    // Create a Molecule instance for testing
    Molecule m("data/geom.dat", 0, 7);

    // Perform tests on the bond function
    double expectedValue = 2.0786985874367461; 
    double tolerance = 1e-5; 

    EXPECT_NEAR(m.bond(0, 1), expectedValue, tolerance);
}

TEST(MoleculeTest, AngleTest) {
    // Create a Molecule instance for testing
    Molecule m("data/geom.dat", 0, 7);

    // Perform tests on the angle function
    double expectedValue = 0.66322511575804732; 
    double tolerance = 1e-5; 
    EXPECT_NEAR(m.angle(0, 1, 2), expectedValue, tolerance);
    // Add more test cases as needed
}

TEST(MoleculeTest, ReadIntegralTest) {
    Molecule m("data/geom.dat", 0, 7);
    m.read_1e_integral("int1e_overlap");
    double expectedValue = 1.0;
    double tolerance = 1e-5; 
    EXPECT_NEAR(m.one_e_overlap_en()[0], expectedValue, tolerance); 
    expectedValue = 0.236703936510848;
    EXPECT_NEAR(m.one_e_overlap_en()[1], expectedValue, tolerance); 
}

// Add more test cases for other functions in a similar manner

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
