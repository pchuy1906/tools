import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read, write
from ase.io.extxyz import read_extxyz
import re, sys
from typing import Dict, List, Any, Optional, Union
import gzip
import pickle
import warnings
import argparse

parser = argparse.ArgumentParser(description="Convert extxyz to pckl.gzip for PACE optimization.")
parser.add_argument('input', help='Input XYZ file')
parser.add_argument("-n", "--name", help="Name dataset for custom weights, if needed.")
args = parser.parse_args()
filename = args.input
configs = read(filename, index=':')
filename = filename.replace(".xyz", ".pckl.gzip") 
all_frames = pd.DataFrame()
reference_energy = 0
for atoms in configs:
    pace_data = {}
    atom_list = atoms.get_chemical_symbols()
    forces = np.array(atoms.get_forces(), dtype = float)
    try:
        stress = atoms.get_stress(voigt=False)
        tstress = True
    except NotImplementedError:
        tstress = False
    energy = atoms.get_potential_energy()
    frame_data = {'energy': [energy],
        'forces': [forces],
        'ase_atoms': [atoms],
        'energy_corrected': [energy - reference_energy],
        'natoms': len(atoms)}
    if (tstress):
        frame_data['stress'] = [stress]
    if (args.name):
        frame_data['name'] = args.name
    dp = pd.DataFrame(frame_data)
    all_frames = pd.concat([all_frames, dp], ignore_index=True)

all_frames.to_pickle(filename, compression='gzip', protocol=4)
#for row in all_frames.itertuples():
#    print(f"natom= {row.natoms}")
#    print(f"name= {row.name}")
