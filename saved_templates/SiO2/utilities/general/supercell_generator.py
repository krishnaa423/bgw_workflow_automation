from pymatgen.core import Lattice, Structure
from pymatgen.io.pwscf import PWInput

# Read CIF file. 
struct = Structure.from_file('SiO2.cif')

# Create pseudo dict, just to generate the pwscf input file. Will not use the dict values. 
pseudo = {
    'Si4+': 'Si-fr.upf',
    'O2-': 'O-fr.upf'
}

# Create pwscf input file for normal cell. 
PWInput(struct, pseudo).write_file('pwscf_normal_cell.in')

# Create supercell CIF file. 
struct.make_supercell([2, 1, 1])
struct.to(filename='SiO2_supecell.cif')
# struct.to(filename='SiO2_supecell.xsf')

# Create pwscf input file for SiO2 supercell. 
PWInput(struct, pseudo).write_file('pwscf_super_cell.in')
