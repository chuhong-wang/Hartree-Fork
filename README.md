# C++ Quantum Chemistry Implementation

This project is a C++ implementation of quantum chemistry calculations from scratch. This package is designed for studying the electronic structure of molecular systems.

## Project Overview

The goal of this project is to provide a C++ version of the core methods in quantum chemistry. This package enables electronic structure calculation of simple systems on a single-CPU. Computing efficiency will be benchmarked against pyscf. 

## Completed Functionality

The current version of the project includes the following features:

- **SCF Convergence:** Iterative self-consistent field (SCF) calculations with convergence criteria.
- **MP2 perturbation energy:** transform AO basis to MO basis and compute MP2 energy. 

## TODO: Work-in-Progress

The project is actively being developed, and the following features are planned for implementation:

- [ ] **Hartree-Fock Calculation:** The Hartree-Fock method is implemented for calculating the electronic ground state of molecular systems.

- [ ] **Geometry Optimization:** Basic geometry optimization functionality to find the equilibrium molecular structure.

- [ ] **Density Functional Theory (DFT):** Implement DFT methods for improved accuracy in electronic structure calculations.

- [ ] **Post-SCF Methods:** Extend the functionality to include post-SCF methods such as correlated wave function methods (e.g., MP2, CCSD).

- [ ] **Parallelization:** Introduce parallel computing capabilities to enhance performance on multi-core systems.

- [ ] **Improved Input/Output:** Enhance input and output handling for improved usability and compatibility.

- [ ] **Documentation:** Provide comprehensive documentation to guide users and developers on utilizing and contributing to the project.

## Getting Started

To build and run the project, follow these steps:

1. Clone the repository:

   ```bash
   git clone git@github.com:chuhong-wang/Hartree-Fork.git
   cd Hartree-Fork 

2. build the project:
    mkdir build
    cd build
    cmake ..
    make

3. Run an example
    TODO 

