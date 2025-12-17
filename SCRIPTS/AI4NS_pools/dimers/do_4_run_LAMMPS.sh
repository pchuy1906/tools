rm -rf dump_xyz_fxyz

pools="$HOME/tools/others/SCRIPTS/AI4NS_pools/"

# file_pairwise is in the $pools
#file_pairwise="lammps_parameters_pairwise_llm105_buck.lmp_zero_ethyl_acetate"
file_pairwise="NEW_PAIR"
#file_pairwise="lammps_parameters_pairwise_llm105_buck.lmp"

# path_lammps_data relatives to the workdir
path_lammps_data="lammps_data_2molecule.lmp"

lammps_inputs=(
        lammps_main.lmp 
    lammps_parameters_covalent_ethyl_acetate.lmp
    lammps_parameters_covalent_llm105.lmp
    lammps_parameters_covalent_pentanol.lmp
    lammps_parameters_covalent_water.lmp
        lammps_forcefield_style.lmp
)


# copy LAMMPS input files:
for file in "${lammps_inputs[@]}"; do
    src="$pools/$file"
    if [[ -f "$src" ]]; then
        cp $src .
    else
        echo "Warning: $src not found"
    fi
done
# change LAMMPS control input file:
sed -i 's|^read_data.*|read_data '"$path_lammps_data"'|'  lammps_main.lmp
# copy pairwise interaction
cp $pools/$file_pairwise lammps_parameters_pairwise_llm105_buck.lmp
# generate lammps_commands.lmp and run LAMMPS
#cp $pools/do_run_LAMMPS.sh .


echo "bond_style     none" >  lammps_commands.lmp
echo "angle_style    none" >> lammps_commands.lmp
echo "dihedral_style none" >> lammps_commands.lmp
echo "improper_style none" >> lammps_commands.lmp
echo >> lammps_commands.lmp

ids=`ls -d VASP_* | sed 's/VASP_//' | sort -n`
for id in $ids ; do

  {
    echo "read_dump dump.scan.lammpstrj                $id x y z"
    echo "fix 1 all nve"
    echo "dump dump1 all custom 1 dump_xyz_fxyz id mol type x y z fx fy fz"
    echo "dump_modify dump1 append yes"
    echo "run 0"
    echo "undump dump1"
    echo "unfix 1"
    echo ""
  }
done >> lammps_commands.lmp

LAMMPSexe="/p/lustre1/pham20/codes/LAMMPS/lammps-4Feb2025/src/lmp_mpi"
${LAMMPSexe} -i lammps_main.lmp
