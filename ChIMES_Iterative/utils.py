import shutil
import re
import os
import subprocess
import socket
from ase.io import read, write
from ase.build import sort
import numpy as np
import glob
from pathlib import Path

pseudo_path = "/usr/workspace/pham20/codes/rzwhippet/VASP_potpaw_PBE"

gpa_2_atm = 9869.2327

machine = socket.gethostname()

if "tuolumne" in machine:
    exe_vasp = "/usr/gapps/qsg/codes/VASP/tuolumne/v5.4.4/vasp_gam"
    exe_lmp = ("/usr/workspace/fuse/lammps/lammps-2Apr2025/"
               "build_jas_tuo_11_18_2025/lmp")
    workflow_type = "FLUX",
    bank = "fuse"
    num_nodes = 5
    cpus_per_task = 96
    max_allow_nodes = 20
elif "dane" in machine:
    exe_vasp = "/usr/gapps/emc-vasp/vasp.6.3.0_vtst/bin/vasp_gam"
    exe_lmp = ("/p/lustre1/pham20/codes/LAMMPS/lammps-4Feb2025/"
               "build_chp_2026_01_25/lmp_mpi")
    workflow_type = "SLURM"
    bank = "pbronze"
    num_nodes = 2
    cpus_per_task = 112
    max_allow_nodes = 8
elif "rzwhippet" in machine:
    exe_vasp = "/usr/gapps/emc-vasp/vasp.6.3.0_vtst/bin/vasp_gam"
    exe_lmp = ("/usr/workspace/pham20/codes/rzwhippet/chimes_calculator_2026_04_17/"
               "etc/lmp/exe/lmp_mpi_chimes")
    workflow_type = "SLURM"
    bank = "ecopper"
    num_nodes = 1
    cpus_per_task = 112
    max_allow_nodes = 20
    stime = "04:00:00"
    queue = "pdebug"
else:
    raise ValueError(f"{machine} is unknown")
nmpi = num_nodes * cpus_per_task



def copy_md_folders(original_dir):
    # 1. Extract the "T" number and generate previous_dir
    match = re.search(r'(T)(\d+)', original_dir)
    if not match:
        raise ValueError(f"Could not find a 'T' number in path: {original_dir}")

    current_t = int(match.group(2))
    previous_dir = original_dir.replace(f"T{current_t}", f"T{current_t - 1}")

    if not os.path.exists(previous_dir):
        raise FileNotFoundError(f"Source directory not found: {previous_dir}")

    # 2. Iterate and copy ONLY data.lammps
    for folder_name in os.listdir(previous_dir):
        prev_folder_path = os.path.join(previous_dir, folder_name)
        source_file = os.path.join(prev_folder_path, "run-1", "data.lammps")

        if os.path.isdir(prev_folder_path) and os.path.exists(source_file):
            # Create the path: original_dir/folder_name/run-1/
            dest_run_dir = os.path.join(original_dir, folder_name, "run-1")

            if not os.path.exists(dest_run_dir):
                os.makedirs(dest_run_dir)

            dest_file = os.path.join(dest_run_dir, "data.lammps")

            if not os.path.exists(dest_file):
                shutil.copy2(source_file, dest_file)
                print(f"Copied data.lammps for {folder_name}")
            else:
                print(f"Skipped {folder_name} (file already exists)")



def write_xyzf(file_VASP_OUTCAR="OUTCAR", file_VASP_POSCAR="POSCAR", cell_type="NON_ORTHO", export_stress=True):

    # parameters
    eV2kcalmol=23.061
    eV2Ha=0.0367493
    Ang2Bohr=1.889725989
    kB2GPa=0.1
    
    #print ("***************************************")
    #print ("reading file POSCAR:", file_VASP_POSCAR)
    f = open(file_VASP_POSCAR, 'rt')
    readPOSCAR = True
    while readPOSCAR:
        for i in range(5):
            line = f.readline()
        line = f.readline().split()
        asym = line
        ntype= len(asym)
        #print (line)
        line = f.readline().split()
        anum = [int(line[i]) for i in range(ntype)]
        #print (line)
        readPOSCAR = False
    f.close()
    
    AtomList = []
    for i in range(ntype):
        ntypei = anum[i]
        for j in range(ntypei):
            AtomList.append(asym[i])
    
    natom = len(AtomList)
    xyz = np.zeros(shape=(natom,3))
    fxyz = np.zeros(shape=(natom,3))
    
    #print ("***************************************")
    #print ("reading file OUTCAR:", file_VASP_OUTCAR)
    f = open(file_VASP_OUTCAR, 'rt')
    while True:
        line = f.readline()
        if line == '': break
    
        keywords = "in kB"
        if keywords in line:
            #print ("*****")
            #print ("read stresses")
            line = line.split()
            #print (line)
            # VASP order is: XX YY ZZ XY YZ ZX
            stresses = [float(line[i])*kB2GPa for i in range(2,8)]
            # ChIMES input wants: XX YY ZZ XY ZX YZ
            tmp = stresses[4]
            stresses[4] = stresses[5]
            stresses[5] = tmp
            str_stresses = ' '.join(map(str, stresses))
            #print (str_stresses)
            done_read_stress = True
    
        keywords = "direct lattice vectors"
        if keywords in line:
            cell_xyz = np.array([])
            #print ("*****")
            #print ("read cell parameter")
            for i in range(3):
                line = f.readline().split()
                #print (line)
                tmp = [float(line[i]) for i in range(0,3)]
                cell_xyz = np.append(cell_xyz, np.array(tmp))
            str_cell_xyz = ' '.join(map(str, cell_xyz))
            #print (str_cell_xyz)
    
        keywords = "POSITION"
        if keywords in line:
            #print ("*****")
            #print ("read atomic coordinates and forces")
            line = f.readline().split()
            for i in range(natom):
                line = f.readline().split()
                #print (line)
                xyz[i,:] = [float(line[i]) for i in range(0,3)]
                fxyz[i,:] = [float(line[i]) for i in range(3,6)]
                ''' VASP has force (eV/A), 
                  while ChIMES used (Ha/Bohr) for input
                  and (kcal/mol/A) for output
                '''
            fxyz = fxyz * eV2Ha/(Ang2Bohr)
    
        keywords = "free  energy"
        if keywords in line:
            #print ("*****")
            #print ("read energy")
            #print (line)
            energy = float(line.split()[4])
            #print (energy)
            energy *= eV2kcalmol
    f.close()
    
    f2 = open("VASP_2_ChIMES.xyzf", "w")
    
    f2.write("%1d\n" %( natom ))
    
    # write the cell parameter
    if (cell_type=="cell_3"):
        for i in range(3):
            j = 4*i
            f2.write("%15.9f" %( cell_xyz[j]))
    elif (cell_type=="NON_ORTHO"):
        f2.write("%s" %( "NON_ORTHO " ))
        for i in range(9):
            f2.write("%15.9f" %( cell_xyz[i]))
    # write stresses
    if export_stress:
        for i in range(6):
            f2.write("%15.9f" %( stresses[i]))
    # write energy
    f2.write("%20.9f" %( energy ))
    
    f2.write("\n")
    for i in range(natom):
        f2.write("%s" %( AtomList[i] ))
        for j in range(3):
            f2.write("%15.9f" %( xyz[i,j]))
        for j in range(3):
            f2.write("%15.9f" %( fxyz[i,j]))
        f2.write("\n")



def write_job_vasp(filename="file_job_vasp"):
    content = f"""#!/bin/sh

#SBATCH -N {num_nodes}
#SBATCH -J dft
#SBATCH -t {stime}
#SBATCH -p {queue}
#SBATCH -A {bank}
#SBATCH --exclusive

nnodes={num_nodes}
nMPI={nmpi}

mv INCAR_PBE INCAR
srun -N $nnodes -n $nMPI {exe_vasp} > out_PBE

mv INCAR_SCAN INCAR
srun -N $nnodes -n $nMPI {exe_vasp} > out_SCAN

touch job_done_vasp"""

    with open(filename, "w") as f:
        f.write(content)


def write_poscar(file_dump="../dump_xyz_vxyz"):
    atoms = read(file_dump, format='lammps-dump-text', index='-1')
    atom_types = atoms.get_array('type')
    mapping = {1: 'C', 2: 'H', 3: 'N', 4: 'O'}
    symbols = [mapping[t] for t in atom_types]
    atoms.set_chemical_symbols(symbols)
    sorted_atoms = sort(atoms)
    write('POSCAR', sorted_atoms, format='vasp')


def get_poscar_elements(file_poscar):
    with open(file_poscar, 'r') as f:
        lines = f.readlines()
        elements = lines[5].split()
        return elements

def write_potcar(file_poscar):
    elements = get_poscar_elements(file_poscar)
    with open("POTCAR", "wb") as pot_out:
        for element in elements:
            path = f"{pseudo_path}/{element}/POTCAR"
            with open(path, "rb") as pot_in:
                pot_out.write(pot_in.read())


def write_kpoints():
    """Writes a Gamma-point KPOINTS file."""
    content = """EA
0
Gamma
1 1 1"""
    with open("KPOINTS", "w") as f:
        f.write(content)



def write_incar(functional_type="PBE"):
    """
    Writes a VASP INCAR file based on the functional type ('PBE' or 'SCAN').
    """
    # Define the differences
    if functional_type.upper() == "SCAN":
        settings = {
            "istart": 2,
            "gga_line": "   METAGGA= SCAN",
            "lwave": ".FALSE.",
            "filename": "INCAR_SCAN"
        }
    else:
        settings = {
            "istart": 0,
            "gga_line": "   GGA    = PE",
            "lwave": ".TRUE.",
            "filename": "INCAR_PBE"
        }

    content = f"""SYSTEM = gas

 Startparameter for this Run:
   NWRITE = 1    LPETIM=F    write-flag & timer
   ISTART = {settings['istart']}    job   : 0-new  1-cont  2-samecut
   NCORE  = 16

 Electronic Relaxation 1
   ENCUT  = 800.0 Set ENCUT to somewhat more than PREC=Accurate value with all elements.
{settings['gga_line']}
   PREC   = Accurate
   EDIFF  = 1E-06  stopping-criterion for ELM
   ISPIN  = 1

 Ionic Relaxation
   NSW    = 0 number of steps for IOM
   IBRION = 2
   ADDGRID= .FALSE.
   ISIF   = 2    stress and relaxation
   LCORR  = .TRUE.    Harris-correction to forces.  Could compare to LCORR=F

 DOS related values:
   ISMEAR = 0;   SIGMA  =    0.04 all-purpose smearing for metals or insulators.

 Electronic Relaxation 2
   NELMIN = 4
   NELM   = 2000
   IALGO  = 38   algorithm
   LREAL  = .FALSE.   NO real space projections / compare to LREAL = A.

   TEBEG  = 0     !Starting temperature
   TEEND  = 0

 Additional Options
  LWAVE   = {settings['lwave']}
  LAECHG  = .FALSE.
  LCHARG  = .FALSE.

LASPH = .TRUE.
ALGO  = Conjugate
NPAR = 9
ISYM=-1"""

    with open(settings['filename'], "w") as f:
        f.write(content)



def write_file_job_lmp(file_job_lmp="job_lmp.sh"):
    content = f"""#!/bin/sh

#SBATCH -N {num_nodes}
#SBATCH -J LMP
#SBATCH -t {stime}
#SBATCH -p {queue}
#SBATCH -A {bank}
#SBATCH --exclusive

nnodes={num_nodes}
nMPI={nmpi}

exe="{exe_lmp}"

srun -N $nnodes -n $nMPI $exe < in.lammps > out.lammps

touch job_done_lmp
"""
    with open(file_job_lmp, "w") as f:
        f.write(content)



def get_total_nodes():
    cmd = "squeue -u $USER"
    output = subprocess.check_output(cmd, shell=True, text=True).strip()
    
    lines = output.split('\n')
    if len(lines) < 2:
        return 0
    
    total_nodes = 0
    for line in lines[1:]:
        parts = line.split()
        if len(parts) >= 7:
            try:
                total_nodes += int(parts[6])
            except ValueError:
                continue
                
    return total_nodes


def write_lammps_input_md(n_md_step=100, temp=3000, n_print=5, file_lmp_input="in.lammps"):
    # Using double braces {{ }} for LAMMPS variables so Python doesn't try to fill them
    content = f"""variable Tini   equal {temp}
variable nprint equal {n_print}

units            real
newton           on

atom_style       charge
atom_modify      map array
atom_modify      sort 0 0.0

read_data        data.lammps
neighbor         1.0 bin
neigh_modify     delay 0 every 1 check no page 100000 one 10000

pair_style       chimesFF
pair_coeff       * * ../../../1-fit/ChIMES_params.txt

velocity         all create ${{Tini}} 12357
fix              1 all nvt temp ${{Tini}} ${{Tini}} 100
thermo_style     custom time step temp etotal pe ke press vol density enthalpy pxx pyy pzz lx ly lz
thermo_modify    format float %20.15g flush yes
thermo           ${{nprint}}
dump             dump_1 all custom ${{nprint}} dump_xyz_vxyz id type x y z vx vy vz
atom_modify      sort 0 0.0

#restart          ${{nprint}} restart.a restart.b

timestep         0.1
run              {n_md_step}
"""
    with open(file_lmp_input, "w") as f:
        f.write(content)


def write_lammps_input_vcopt(p_gpa=0, file_lmp_input="in.lammps"):
    p_atm = p_gpa * gpa_2_atm
    content = f"""variable p0 equal {p_atm}

units            real
newton           on

atom_style       charge
atom_modify      map array
atom_modify      sort 0 0.0

read_data        data.lammps
neighbor         1.0 bin
neigh_modify     delay 0 every 1 check no page 100000 one 10000

pair_style       chimesFF
pair_coeff       * * ../../../1-fit/ChIMES_params.txt

thermo_style     custom time step temp etotal pe ke press vol density enthalpy pxx pyy pzz lx ly lz
thermo_modify    format float %20.15g flush yes
thermo           1
atom_modify      sort 0 0.0

fix              1 all box/relax aniso ${{p0}}
minimize         1.0e-6 1.0e-6 100000 100000

write_dump       all custom dump_xyz_vxyz id type x y z vx vy vz
"""
    with open(file_lmp_input, "w") as f:
        f.write(content)



def get_p_t(folder):
    if "ambient" in folder:
        p, t = 0.0, 300.0
    elif "condensed_new" in folder:
        p, t = 0.0, 3000.0
    elif "condensed_TATB_extreme" in folder:
        t = folder.split('_')[3].replace('K', '')
        p = 0.0
    elif "Temp" in folder:
        t = folder.split('_')[-1]
        p = 0.0
    elif "vcopt" in folder:
        p = folder.split('_')[-2]
        t = 0.0
    else:
        raise ValueError(f"Unexpected folder format: '{folder}'")
    return float(p), float(t)

def clean_vasp():
    list_files = ["CHG", "CHGCAR", "DOSCAR", "EIGENVAL", "IBZKPT", "OSZICAR", "PCDAT", "vasprun.xml", "WAVECAR", "XDATCAR"]
    [Path(f).unlink(missing_ok=True) for f in list_files]
    for f in glob.glob("slurm*"):
        os.remove(f)

