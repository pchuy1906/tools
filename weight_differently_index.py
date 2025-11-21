import os
import sys

import numpy as np

# read input 
import argparse
parser = argparse.ArgumentParser(description='generate weights')
# Arguments supported by the code.
parser.add_argument("--nAtomType",           type=int,   default=4)
parser.add_argument("--nCondensed",          type=int,   default=10)
parser.add_argument("--wE",                  type=float, default=200.0)
parser.add_argument("--wS",                  type=float, default=500.0)

parser.add_argument("--wE_more",             type=float, default=5.0)
parser.add_argument('--file_id_more_weight',             default="aaa.dat")

parser.add_argument("--more_force_H",                    default=False, action="store_true", help="weight more on forces of hydrogen atoms")
parser.add_argument('--file_xyzf',                       default="aaa.xyzf")
parser.add_argument("--w_force_H_more",      type=float, default=5.0)

args    = parser.parse_args()
nAtomType           = args.nAtomType
nCondensed          = args.nCondensed
wE                  = args.wE
wS                  = args.wS
wE_more             = args.wE_more
file_id_more_weight = args.file_id_more_weight
more_force_H        = args.more_force_H
file_xyzf           = args.file_xyzf
w_force_H_more      = args.w_force_H_more

if not os.path.isfile(file_id_more_weight):
    print(f"ERROR: {file_id_more_weight} does not exist.")
    sys.exit(1)
else:
    print(f"configurations with IDs in file {file_id_more_weight} is weighted more with a factor {wE_more}")
ids_more_weight = np.loadtxt(file_id_more_weight, dtype=int)

f2 = open('label.txt', "w")
f3 = open('new_weight.dat', "w")

def read_xyz_trajectory_atomic_symbols(filename):
    """
    Reads an XYZ trajectory file and extracts atomic symbols for each configuration.

    Args:
        filename (str): The path to the XYZ trajectory file.

    Returns:
        list: A list of lists, where each inner list contains the atomic symbols
              for a single configuration (frame) in the trajectory.
    """
    all_atomic_symbols = []
    with open(filename, 'r') as f:
        while True:
            # Read the number of atoms for the current frame
            num_atoms_line = f.readline()
            if not num_atoms_line:  # End of file
                break
            num_atoms = int(num_atoms_line.strip())

            # Read the comment line (ignored)
            f.readline()

            # Read atomic symbols for the current frame
            current_frame_symbols = []
            for _ in range(num_atoms):
                line = f.readline()
                if not line: # Handle unexpected end of file
                    break 
                parts = line.split()
                if parts:
                    current_frame_symbols.append(parts[0])
            
            if current_frame_symbols: # Only add if symbols were found for the frame
                all_atomic_symbols.append(current_frame_symbols)

    return all_atomic_symbols

if more_force_H:
    if not os.path.isfile(file_xyzf):
        print(f"ERROR: {file_xyzf} does not exist.")
        sys.exit(1)
    else:
        print(f"force on hydrogen atoms will be weighted more, with w_force_H_more = {w_force_H_more} ")
atomic_symbols_trajectory = read_xyz_trajectory_atomic_symbols(file_xyzf)

def gen_molecule_name(AType, ANumber):
    mole_name = ""
    nAtomType = len(AType)
    for i in range(nAtomType):
        mole_name = mole_name + AType[i] + str(ANumber[i])
    return mole_name

tmp_file = "frames.all.log"

f  = open(tmp_file ,"rt")
while True:

    tmp  = f.readline()
    # if EOF, stop 
    line = tmp.strip()
    if line == '': break
    # otherwise, continue reading file
    # The first line look like this: Processing frame 0 on rank 0
    tmp  = tmp.split()
    id_frame = int(tmp[2])

    tmp  = f.readline().split()
    AType = []
    ANumber = []
    for j in range(nAtomType):
        AType.append(tmp[2*j])
        ANumber.append(int(tmp[2*j+1]))
    mole_name = gen_molecule_name(AType, ANumber)

    tmp  = f.readline()

    natom = sum(ANumber)
    if more_force_H:
        if natom != len(atomic_symbols_trajectory[id_frame]):
            sys.exit(1)

    for j in range(natom):
        f2.write("force_"+atomic_symbols_trajectory[id_frame][j]+"_"+mole_name+"\n")
        f2.write("force_"+atomic_symbols_trajectory[id_frame][j]+"_"+mole_name+"\n")
        f2.write("force_"+atomic_symbols_trajectory[id_frame][j]+"_"+mole_name+"\n")
        if more_force_H:
            if atomic_symbols_trajectory[id_frame][j] == "H":
                f3.write('%12.2f \n' %(w_force_H_more) )
                f3.write('%12.2f \n' %(w_force_H_more) )
                f3.write('%12.2f \n' %(w_force_H_more) )
            else:
                f3.write("1.0\n")
                f3.write("1.0\n")
                f3.write("1.0\n")

    if (id_frame < nCondensed):
        f2.write("dia_stress_" +mole_name+"\n")
        f2.write("stress_xy_"  +mole_name+"\n")
        f2.write("stress_xz_"  +mole_name+"\n")
        f2.write("stress_xy_"  +mole_name+"\n")
        f2.write("dia_stress_" +mole_name+"\n")
        f2.write("stress_yz_"  +mole_name+"\n")
        f2.write("stress_xz_"  +mole_name+"\n")
        f2.write("stress_yz_"  +mole_name+"\n")
        f2.write("dia_stress_" +mole_name+"\n")

        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )
        f3.write('%15.6f \n' %(wS) )


    fE_extra = 1.0
    if (id_frame in ids_more_weight):
        fE_extra = wE_more

    f2.write("energy_"+mole_name+"\n")
    f2.write("energy_"+mole_name+"\n")
    f2.write("energy_"+mole_name+"\n")
    f3.write('%15.6f \n' %(fE_extra * wE/float(natom)))
    f3.write('%15.6f \n' %(fE_extra * wE/float(natom)))
    f3.write('%15.6f \n' %(fE_extra * wE/float(natom)))

