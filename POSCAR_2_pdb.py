from ase.io import read, write
atoms = read('POSCAR')
write('structure.pdb', atoms)
