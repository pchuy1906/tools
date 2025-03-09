echo
echo "usage: \"do_script.sh  file.xyz  option-cell(1-cell_3, 2-cell_9, 3-non-orthor)\""
echo

if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

file_xyz="$1"
option_cell="$2"

# CC : 1.90
rcut=1.90
~/tools/others/do_bond_counts_each.sh ${file_xyz} ${option_cell} $rcut    C C &> OUTPUT_bond

# CN : 1.80
rcut=1.80
~/tools/others/do_bond_counts_each.sh ${file_xyz} ${option_cell} $rcut    C N &> OUTPUT_bond

# CO : 1.80
rcut=1.80
~/tools/others/do_bond_counts_each.sh ${file_xyz} ${option_cell} $rcut    C O &> OUTPUT_bond

# NN : 1.75
rcut=1.75
~/tools/others/do_bond_counts_each.sh ${file_xyz} ${option_cell} $rcut    N N &> OUTPUT_bond

# NO : 1.65
rcut=1.65
~/tools/others/do_bond_counts_each.sh ${file_xyz} ${option_cell} $rcut    N O &> OUTPUT_bond

# OO : 1.70
rcut=1.70
~/tools/others/do_bond_counts_each.sh ${file_xyz} ${option_cell} $rcut    O O &> OUTPUT_bond


