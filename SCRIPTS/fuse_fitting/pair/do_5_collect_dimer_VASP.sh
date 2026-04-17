cat > ddd.py << EOF
import sys
import numpy as np

file_VASP_energy = sys.argv[1]
energies = np.loadtxt(file_VASP_energy)

each_energy = np.array([energies[-2], energies[-1]])
each_n_molecule = np.loadtxt("decomposition.dat")

total_energy = np.dot(each_energy.flatten(), each_n_molecule)

eV_2_kcalmol = 23.0609
interaction_energy = (energies-total_energy) * eV_2_kcalmol
np.savetxt("EINT.txt", interaction_energy[:-2], fmt="%.6f")
EOF




cwd=`pwd`
for fold in $(ls -1vd VaspOracle/dimer* ) ; do
    cd $fold
        path=`pwd`
        fold_name=`echo $path | sed 's|/| |g' | awk '{print $NF}'`

        file="VASP_ENERGY.dat"
        rm -rf $file
        for fileOUTCAR in $(ls -1vd */OUTCAR) ; do
            grep -a "free  energy" $fileOUTCAR | awk '{print $5}' >> $file
        done

        echo "1 1" > decomposition.dat

        python $cwd/ddd.py $file

    cd $cwd
done




