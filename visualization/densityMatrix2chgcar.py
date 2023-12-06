import numpy as np
from pyscf import gto, scf
from pyscf.tools import chgcar
# cell = gto.M(atom='H 0 0 0; H 0 0 1', a=numpy.eye(3)*3)
# mf = scf.RHF(cell).run()


class Tools:
    def __init__(self):
        self.atom_charges = None  # type: np.ndarray
        self.atom_coords = None  # type: np.ndarray
        self.natm = None  # type: int
        self.mol = None # type: gto.mole 
        self.dm = None # type: np.ndarray 
        self.charge = 0
        
    def mol_instance_from_dat_file(self, file_path: str, basis:str):
        with open(file_path, "r") as f:
            dat = np.array([line.split() for line in f.readlines()][1:])
            self.atom_charges = np.array(dat[:, 0], dtype=float).astype(int)
            self.atom_coords = np.array(dat[:, 1:4], dtype=float)
            self.natm = self.atom_charges.shape[0]
        mol = gto.Mole()
        mol.unit = "Bohr"
        mol.atom = "\n".join([("{:3d} " + " ".join(["{:25.18f}"] * 3)).format(chg, *coord) for chg, coord in zip(self.atom_charges, self.atom_coords)])
        mol.basis = basis
        mol.charge = self.charge
        # mol.spin = 0
        mol.multiplicity = 1
        self.mol = mol.build()
    
    def read_density_matrix(self, file_path:str):
        with open(file_path, "r") as f:
            dat = np.array([line.split() for line in f.readlines()])
            self.dm = np.array(dat[:,:], dtype=float)
            
    
    def dm2ed(self, outfile:str):
        print("shape of dm", self.dm.shape)
        chgcar.density(self.mol, outfile, self.dm)
        
for i in range(10):
    t=Tools()
    t.mol_instance_from_dat_file("data/geom.dat", "STO-3G")
    t.read_density_matrix("data/dm_"+str(i)+".dat")
    t.dm2ed(str(i)+".CHGCAR")
