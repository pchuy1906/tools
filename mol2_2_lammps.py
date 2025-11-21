import parmed
print(dir(parmed.formats))

print(parmed.__version__)
structure = parmed.amber.AmberParm('structure.prmtop', 'structure.inpcrd')
structure.save('structure.lammps', format='lammps')
