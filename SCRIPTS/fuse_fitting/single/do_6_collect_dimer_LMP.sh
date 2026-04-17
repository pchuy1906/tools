cat > eee.py << EOF
import sys
import numpy as np
from orchestrator.target_property.analysis import AnalyzeLammpsLog

file_lammps_log = sys.argv[1]
lammps_log = AnalyzeLammpsLog(file_lammps_log)
energies = lammps_log.get('PotEng')

file_lammps_log = sys.argv[2]
lammps_log = AnalyzeLammpsLog(file_lammps_log)
energy_1 = lammps_log.get('PotEng')

file_lammps_log = sys.argv[3]
lammps_log = AnalyzeLammpsLog(file_lammps_log)
energy_2 = lammps_log.get('PotEng')

each_energy = np.array([energy_1, energy_2])
each_n_molecule = np.loadtxt("decomposition.dat")

total_energy = np.dot(each_energy.flatten(), each_n_molecule)

interaction_energy = energies-total_energy
np.savetxt("EINT.txt", interaction_energy, fmt="%.6f")
EOF




cwd=`pwd`
for fold in $(ls -1vd LAMMPSSimulator/dimer_* ) ; do
    cd $fold
        pwd
        echo "1 1" > decomposition.dat
        python $cwd/ccc.py  00000/lammps.out 00001/lammps.out 00002/lammps.out

    cd $cwd
done




