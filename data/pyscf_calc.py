import numpy as np
from pyscf import gto, scf
from pyscf.tools import chgcar

def calc_idx(i,j):
    if i>j: return int((i+1)*i/2+j)
    else: return int((j+1)*j/2+i)

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
            lines = f.readlines()
            line0 = lines[0]
            dat = np.array([line.split() for line in lines][1:])
            self.atom_charges = np.array(dat[:, 0], dtype=float).astype(int)
            self.atom_coords = np.array(dat[:, 1:4], dtype=float)
            # self.atom_coords*= 1.88973 # angstrom to Bohr
            self.natm = self.atom_charges.shape[0]
        # with open(file_path, "w") as f:
        #     f.write(line0)
        #     lines = np.concatenate((self.atom_charges.reshape(-1, 1), self.atom_coords), axis=1)
        #     f.write('\n'.join(' '.join(map(str, line)) for line in lines))

        mol = gto.Mole()
        mol.unit = "Bohr"
        mol.atom = "\n".join([("{:3d} " + " ".join(["{:25.18f}"] * 3)).format(chg, *coord) for chg, coord in zip(self.atom_charges, self.atom_coords)])
        mol.basis = basis
        mol.charge = self.charge
        # mol.spin = 0
        mol.multiplicity = 1
        self.mol = mol.build()

t = Tools()
t.mol_instance_from_dat_file('geom.dat', basis='STO-3G')

intor_dict = {'int1e_ovlp':"s.dat", 'int1e_kin':"t.dat", 'int1e_nuc':"v.dat", "int2e":"eri.dat"}
mol = t.mol

for inte_key in intor_dict:
    inte = mol.intor(inte_key)
    output = ''
    for i in range(len(inte)):
        for j in range(len(inte[0])):
            if i>=j:
                output+= str(i+1) + " " + str(j+1) + " " + str(inte[i][j]) + '\n'
    
    with open(intor_dict[inte_key], 'w') as f:
        f.write(output)

output = ''
inte_key = 'int2e'
inte = mol.intor(inte_key)
print(inte.shape, inte[0][0][0][0])
for i in range(len(inte)):
        for j in range(len(inte[0])):
            if i>=j:
                for k in range(len(inte[0][0])):
                    for l in range(len(inte[0][0][0])):
                        output+= str(i+1) + " " + str(j+1) + " " + str(k+1) + " " + str(l+1) + " " + str(inte[i][j][k][l]) + '\n'
with open(intor_dict[inte_key], 'w') as f:
    f.write(output)

from pyscf.gto.mole import classical_coulomb_energy
f = open('enuc.dat', 'w')
f.write(str(classical_coulomb_energy(mol)))
f.close()