import argparse
parser = argparse.ArgumentParser(description='traj.ase to traj.xyz')
parser.add_argument('--cell_9', nargs='+', type=float)
parser.add_argument('--traj_ase', default="traj.ase")
args = parser.parse_args()
traj_ase = args.traj_ase
cell_9   = args.cell_9


import ase
from ase.io import read, write
from ase.io.vasp import write_vasp_xdatcar
from ase.io.trajectory import Trajectory

bec = Trajectory(traj_ase)
bec_new = ase.io.write("tmp.xyz", bec, format="xyz")


f  = open("tmp.xyz" ,"rt")
f2 = open("result_traj.xyz", "w")

while True:
    tmp  = f.readline()
    f2.write("%s" %( tmp ))

    line = tmp.strip()
    if line == '': break

    natom = int(tmp.split()[0])

    tmp  = f.readline()
    f2.write("%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n" %( cell_9[0], cell_9[1], cell_9[2], cell_9[3], cell_9[4], cell_9[5], cell_9[6], cell_9[7], cell_9[8] ))

    for i in range(natom):
        tmp  = f.readline()
        f2.write("%s" %( tmp ))


