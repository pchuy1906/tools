from ase.io import read, write

atoms = read("OUTCAR", format='vasp-out')

#print(dir(atoms))
#print (atoms.calc)
#print(dir(atoms.calc))
#print(atoms.calc.results['free_energy'])

stress = -1.0*atoms.get_stress(voigt=False)
atoms.calc.results['stress'] = stress

free_energy = atoms.calc.results['free_energy']
atoms.calc.results['energy'] = free_energy
del atoms.calc.results['free_energy']

write("OUTCAR.extxyz", atoms, format='extxyz')
