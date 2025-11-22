#!/usr/bin/env python3
"""
Script to set up VASP calculations from LAMMPS dump and data files.
"""

import os
import sys
import argparse
import logging
import numpy as np
from typing import Sequence, Any

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Setup VASP calculations from LAMMPS dump and data files.'
    )
    parser.add_argument(
        '--file_LAMMPS_dump',
        default='dump.lammps',
        help='Path to the LAMMPS dump file'
    )
    parser.add_argument(
        '--file_LAMMPS_data',
        default='data.lammps',
        help='Path to the LAMMPS data file'
    )
    return parser.parse_args()

## Arguments supported by the code.
#parser.add_argument("--file_dump",   default='file.dump', help='file with format dump lammps')
#parser.add_argument("--file_energy", default='file.ener', help='file with format log lammps')
#parser.add_argument("--export_quan", default='xyzfes',    help='export xyzfe/xyzfes')
#parser.add_argument("--cell_type",   default='NON_ORTHO', help='cell_type (cell_3/cell_9/NON_ORTHO)')
#parser.add_argument('--atom_types',  nargs='+')
#parser.add_argument("--lmp_unit",    default='real',    help='real/metal')
#parser.add_argument("--read_charge", default=False, action="store_true", help="read point charge")
#
#args    = parser.parse_args()
#file_dump     = args.file_dump
#file_energy   = args.file_energy
#atom_types    = args.atom_types
#export_quan   = args.export_quan
#cell_type     = args.cell_type
#lmp_unit      = args.lmp_unit
#read_charge   = args.read_charge
#
#atm_2_GPa = 0.000101325
#bar_2_GPa = 0.0001
#
#A_2_Bohr = 1.889725989
#kcalmol_2_Ha = 0.00159362
#kcalmol_A_2_Ha_Bohr = kcalmol_2_Ha / A_2_Bohr
#
#eV_2_Ha = 0.0367502 
#eV_A_2_Ha_Bohr = eV_2_Ha / A_2_Bohr
#
#eV_2_kcalmol = 23.0609 
#
#
#f  = open(file_dump ,"rt")
#nconf = 1
#
#fname = "charges.dat"
#f2 = open( fname, "w")
#
#if (export_quan=="xyzfes"):
#    fname = "LAMMPS.xyzfes"
#elif (export_quan=="xyzfe"):
#    fname = "LAMMPS.xyzfe"
#else:
#    print ("unknown option")
#    exit()
#f3 = open( fname, "w")
#
#f0  = open(file_energy ,"rt")
#tmp = f0.readline().split()
#tmp = np.array(tmp)
#ipe  = np.where(tmp=='PotEng')[0][0]
#ixx  = np.where(tmp=='Pxx')[0][0]
#iyy  = np.where(tmp=='Pyy')[0][0]
#izz  = np.where(tmp=='Pzz')[0][0]
#ixy  = np.where(tmp=='Pxy')[0][0]
#ixz  = np.where(tmp=='Pxz')[0][0]
#iyz  = np.where(tmp=='Pyz')[0][0]
#
#while True:
#
#    tmp  = f.readline()
#    line = tmp.strip()
#    if line == '': break
#    tmp  = f.readline()
#    tmp  = f.readline()
#    tmp  = f.readline()
#    natom = int(tmp)
#
#    tmp = f.readline()
#    if "xy" in tmp:
#        #print ("NON_ORTHO")
#
#        tmp = f.readline().split()
#        [xlo_bound, xhi_bound, xy] = [float(x) for x in tmp]
#
#        tmp = f.readline().split()
#        [ylo_bound, yhi_bound, xz] = [float(x) for x in tmp]
#
#        tmp = f.readline().split()
#        [zlo_bound, zhi_bound, yz] = [float(x) for x in tmp]
#
#        zlo = zlo_bound
#        zhi = zhi_bound
#        ylo = ylo_bound - min(0.0,yz)
#        yhi = yhi_bound - max(0.0,yz)
#        xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
#        xhi = xhi_bound - max(0.0,xy,xz,xy+xz)
#
#        lx = xhi-xlo
#        ly = yhi-ylo
#        lz = zhi-zlo
#    else:
#        tmp = f.readline().split()
#        #print (tmp)
#        lx = float(tmp[1])-float(tmp[0])
#        tmp = f.readline().split()
#        #print (tmp)
#        ly = float(tmp[1])-float(tmp[0])
#        tmp = f.readline().split()
#        #print (tmp)
#        lz = float(tmp[1])-float(tmp[0])
#        #print (lx,ly,lz)
#        xy = 0.0
#        xz = 0.0
#        yz = 0.0
#
#    cell_9 = [lx,0,0, xy,ly,0, xz,yz,lz]
#
#    tmp = f.readline().split()
#    tmp = np.array(tmp)
#    if (read_charge):
#        idq  = np.where(tmp=='q')[0][0]
#    idfx = np.where(tmp=='fx')[0][0]
#    idfy = np.where(tmp=='fy')[0][0]
#    idfz = np.where(tmp=='fz')[0][0]
#    idx  = np.where(tmp=='x')[0][0]
#    idy  = np.where(tmp=='y')[0][0]
#    idz  = np.where(tmp=='z')[0][0]
#    idi  = np.where(tmp=='type')[0][0]
#
#    #print (idq, idfx, idfy, idfz, idx, idy, idz, idi)
#    q  = np.zeros(shape=(natom))
#    x  = np.zeros(shape=(natom))
#    y  = np.zeros(shape=(natom))
#    z  = np.zeros(shape=(natom))
#    fx = np.zeros(shape=(natom))
#    fy = np.zeros(shape=(natom))
#    fz = np.zeros(shape=(natom))
#    strs = ["XXXX" for x in range(natom)]
#    atype = np.array(strs)
#    for i in range(natom):
#        tmp = f.readline().split()
#        idt = int(tmp[0])
#        if (read_charge):
#            q[idt-1]  = float(tmp[idq-2])
#        x[idt-1]  = float(tmp[idx-2])
#        y[idt-1]  = float(tmp[idy-2])
#        z[idt-1]  = float(tmp[idz-2])
#
#        f_convert = 1.0
#        p_convert = 1.0
#        e_convert = 1.0
#        if (lmp_unit == "real"):
#            f_convert = kcalmol_A_2_Ha_Bohr
#            p_convert = atm_2_GPa
#            e_convert = 1.0
#        elif (lmp_unit == "metal"):
#            f_convert = eV_A_2_Ha_Bohr
#            p_convert = bar_2_GPa
#            e_convert = eV_2_kcalmol
#
#
#        fx[idt-1] = float(tmp[idfx-2]) * f_convert
#        fy[idt-1] = float(tmp[idfy-2]) * f_convert
#        fz[idt-1] = float(tmp[idfz-2]) * f_convert
#        iid = int(tmp[idi-2])
#        atype[idt-1] = atom_types[iid-1]
#
#    TPM  = f0.readline().split()
#    energy = float(TPM[ipe]) * e_convert
#
#    Pxx = float(TPM[ixx]) * p_convert 
#    Pyy = float(TPM[iyy]) * p_convert 
#    Pzz = float(TPM[izz]) * p_convert 
#    Pxy = float(TPM[ixy]) * p_convert 
#    Pxz = float(TPM[ixz]) * p_convert 
#    Pyz = float(TPM[iyz]) * p_convert 
#
#    f3.write("%d\n" %(natom))
#    if (export_quan=="xyzfes"):
#        if (cell_type=="NON_ORTHO"):
#            f3.write("%s" %("NON_ORTHO" ))
#            #f3.write("%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f" %( lx,0,0, 0,ly,0, 0,0,lz ))
#            for i in range(9):
#                f3.write("%15.6f " %( cell_9[i]))
#            f3.write("%20.6f %20.6f %20.6f %20.6f %20.6f %20.6f" %( Pxx, Pyy, Pzz, Pxy, Pxz, Pyz ))
#            f3.write("%20.6f" %( energy ))
#            f3.write("\n")
#        elif (cell_type=="cell_3"):
#            f3.write("%15.6f %15.6f %15.6f" %( lx,ly,lz ))
#            f3.write("%20.6f %20.6f %20.6f %20.6f %20.6f %20.6f" %( Pxx, Pyy, Pzz, Pxy, Pxz, Pyz ))
#            f3.write("%20.6f" %( energy ))
#            f3.write("\n")
#        else:
#            print ("unknown option")
#            exit()
#    elif (export_quan=="xyzfe"):
#        f3.write("%15.6f %15.6f %15.6f" %( lx,ly,lz ))
#        f3.write("%20.6f" %( energy ))
#        f3.write("\n")
#    else:
#        print ("unknown option")
#        exit()
#
#    if (read_charge):
#        f2.write("%d " %(nconf))
#        for i in range(natom):
#            f2.write("%6.3f " %(q[i]))
#        f2.write("\n")
#
#    for i in range(natom):
#        f3.write("%s" %(atype[i]))
#        f3.write("%15.9f" %(x[i]))
#        f3.write("%15.9f" %(y[i]))
#        f3.write("%15.9f" %(z[i]))
#        f3.write("%15.9f" %(fx[i]))
#        f3.write("%15.9f" %(fy[i]))
#        f3.write("%15.9f" %(fz[i]))
#        f3.write("\n")
#
#
#    nconf += 1
#
#
#f2.close()
#
#

def extract_masses(filename):
    masses = []
    f = open(filename, 'rt')
    while True:
        line = f.readline()
        if line == '': break
        if "atom types" in line:
            ntype = int(line.split()[0])
        if "Masses" in line:
            line = f.readline() # This line is empty
            for i in range(ntype):
                mass = f.readline().split()[1]
                masses.append(float(mass))
    f.close()
    return np.array(masses)

def read_dump(filename, element_symbols):
    with open(filename, "rt") as f:
        nconf = 0
        while True:
            line = f.readline()
            if not line:
                break  # End of file

            # Skip three lines (assuming header or irrelevant info)
            for _ in range(3):
                tmp = f.readline()

            natom_line = tmp.strip()
            try:
                natom = int(natom_line.split()[0])
            except (IndexError, ValueError) as e:
                logging.error(f"Failed to parse atom count: {natom_line} ({e})")
                print(f"Error: Failed to parse atom count: {natom_line}")
                sys.exit(1)

            #logging.info(f"Config {nconf}: Atom count = {natom}")

            # Read cell lines
            cell_info_line = f.readline()
            line_cell_x = f.readline().split()
            line_cell_y = f.readline().split()
            line_cell_z = f.readline().split()

            if "xy" in cell_info_line:
                cell_type = "cell_9"
            else:
                cell_type = "cell_3"

            cell9 = read_cell(cell_type, line_cell_x, line_cell_y, line_cell_z)
            #logging.info(f"Cell parameters: {cell9}")

            # Read atom data header
            atom_header = f.readline().split()
            atom_header = np.array(atom_header) 
            try:
                idx = np.where(atom_header == 'x')[0][0]
                idy = np.where(atom_header == 'y')[0][0]
                idz = np.where(atom_header == 'z')[0][0]
            except IndexError:
                idx = np.where(atom_header == 'xu')[0][0]
                idy = np.where(atom_header == 'yu')[0][0]
                idz = np.where(atom_header == 'zu')[0][0]
            idtype  = np.where(atom_header=='type')[0][0]
            idmol   = np.where(atom_header=='mol')[0][0]

            # Read atom data lines
            atom_data = []
            x  = np.zeros(shape=(natom))
            y  = np.zeros(shape=(natom))
            z  = np.zeros(shape=(natom))
            atype  = np.zeros(shape=(natom), dtype=int)
            amol  = np.zeros(shape=(natom), dtype=int)

            for i in range(natom):
                atom_line = f.readline().split()
                idt = int(atom_line[0])
                x[idt-1]  = float(atom_line[idx-2])
                y[idt-1]  = float(atom_line[idy-2])
                z[idt-1]  = float(atom_line[idz-2])
                atype[idt-1]  = float(atom_line[idtype-2])
                amol[idt-1]  = float(atom_line[idmol-2])

            atomList = element_symbols[atype-1]
            nconf += 1

def main():
    setup_logging()
    args = parse_arguments()

    lammps_dump_path = args.file_LAMMPS_dump
    lammps_data_path = args.file_LAMMPS_data

    logging.info('Starting collecting XYZ from LAMMPS dump format.')

    # Check if files exist
    if not os.path.isfile(lammps_data_path):
        logging.error(f'LAMMPS data file not found: {lammps_data_path}')
        return
    if not os.path.isfile(lammps_dump_path):
        logging.error(f'LAMMPS dump file not found: {lammps_dump_path}')
        return

    # Your conversion logic here
    logging.info(f'Reading masses from LAMMPS data file: {lammps_data_path}')
    masses = extract_masses(lammps_data_path)
    #print (masses)
    #print(np.array2string(masses, separator=', '))
    print(', '.join([f'{m:.3f}' for m in masses]))
    logging.info(f'Matching these masses with element symbols')
    element_symbols = match_mass_to_element_ase(masses)
    print (element_symbols)

    logging.info(f'Reading LAMMPS dump file: {lammps_dump_path}')
    read_dump(lammps_dump_path, element_symbols)


if __name__ == '__main__':
    main()

