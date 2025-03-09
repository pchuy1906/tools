from ase import units
from ase import Atoms
from ase.optimize import BFGS
#from ase.filters import Filter
from ase.vibrations import Vibrations
from ase.io.vasp import read_vasp
from ase.calculators.chimes_ase import ChIMES 

import argparse
parser = argparse.ArgumentParser(description='calculator for ChIMES/ANI1xnr in ase')
# Arguments supported by the code.
parser.add_argument("--file_ChIMES_param", default='',       help='file ChIMES parameter')
parser.add_argument("--file_POSCAR",       default='POSCAR', help='file POSCAR')
parser.add_argument("--calc_type",         default='sp', help='calc type: sp(single-point), opt(optimize)')

args        = parser.parse_args()
file_POSCAR        = args.file_POSCAR
file_ChIMES_param  = args.file_ChIMES_param
calc_type          = args.calc_type

bar_2_GPa = 0.0001
eV_2_kcalmol = 23.0609

eV_2_Ha = 0.0367502
Ang_2_Bohr = 1.8897259886

atoms = read_vasp(file=file_POSCAR)

print (len(file_ChIMES_param))
if (len(file_ChIMES_param)==0):
    import torchani
    calc = torchani.models.ANI1xnr().ase()
else:
    calc = ChIMES(file_ChIMES_param)

atoms.set_calculator(calc)

if (calc_type=="sp"):
    e = atoms.get_potential_energy()
    # ChIMES uses kcal/mol
    e_kcalmol = e * eV_2_kcalmol
    print('Energy', e_kcalmol)
    f = atoms.get_forces()
    f_Ha_Bohr = f * eV_2_Ha / Ang_2_Bohr
    # ChIMES uses Ha/Bohr
    print('Forces')
    print(f_Ha_Bohr)
    
    stress = atoms.get_stress(include_ideal_gas=False)
    # ChIMES uses GPa
    stress_GPa = -stress / units.bar * bar_2_GPa
    print (stress_GPa)
    # Voigt order (xx, yy, zz, yz, xz, xy)
    sig_xx = stress_GPa[0]
    sig_yy = stress_GPa[1]
    sig_zz = stress_GPa[2]
    sig_yz = stress_GPa[3]
    sig_xz = stress_GPa[4]
    sig_xy = stress_GPa[5]
    
    xyz = atoms.get_positions()
    print (xyz)
    symbols =  atoms.get_chemical_symbols()
    print (symbols)
    
    cell =  atoms.get_cell()
    print (cell)
    
    f2 = open("input.xyzfes", "w")
    
    natom = len(symbols)
    f2.write("%d\n" %( natom ))
    
    f2.write("%s" %( "NON_ORTHO " ))
    for i in range(3):
        for j in range(3):
            f2.write("%15.9f" %( cell[i][j] ))
    
    f2.write("%15.9f" %( sig_xx ))
    f2.write("%15.9f" %( sig_yy ))
    f2.write("%15.9f" %( sig_zz ))
    f2.write("%15.9f" %( sig_xy ))
    f2.write("%15.9f" %( sig_yz ))
    f2.write("%15.9f" %( sig_xz ))
    
    f2.write("%22.9f" %( e_kcalmol ))
    f2.write("\n")
    
    natom = len(symbols)
    for i in range(natom):
        f2.write("%4s" %( symbols[i] ))
        for j in range(3):
            f2.write("%15.9f" %( xyz[i,j] ))
        fi = f_Ha_Bohr[i,:]
        for j in range(3):
            f2.write("%15.9f" %( fi[j] ))
        f2.write("\n")


if (calc_type=="opt"):
    dyn = BFGS(atoms, trajectory='opt_traj.ase')
    dyn.run(fmax=0.001)
    e = atoms.get_potential_energy()
    print ('final energy (eV)', e)
    # ChIMES uses kcal/mol
    e_kcalmol = e * eV_2_kcalmol
    print('final energy (kcal/mol)', e_kcalmol)


if (calc_type=="vcopt"):
    print ("problem with Filter class")
    #ucf = UnitCellFilter(atoms)
    #dyn = BFGS(ucf)
    #dyn.run(fmax=0.001)
    #e = atoms.get_potential_energy()
    #print ('final energy (eV)', e)
    ## ChIMES uses kcal/mol
    #e_kcalmol = e * eV_2_kcalmol
    #print('final energy (kcal/mol)', e_kcalmol)



