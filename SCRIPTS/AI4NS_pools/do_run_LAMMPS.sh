fileinputdump="../all_dump"

ids=`grep -A1 TIMESTEP $fileinputdump | awk '/^[0-9]+$/'`

rm -rf dump_xyz_fxyz

echo "bond_style     none" >  lammps_commands.lmp
echo "angle_style    none" >> lammps_commands.lmp
echo "dihedral_style none" >> lammps_commands.lmp
echo "improper_style none" >> lammps_commands.lmp
echo >> lammps_commands.lmp

for id in $ids ; do

  {
    echo "read_dump    $fileinputdump                $id x y z"
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




rm -rf VASP_2_ChIMES.xyzf
python ~/tools/others/AI4NS_from_dump.py \
  --file_LAMMPS_dump dump_xyz_fxyz  \
  --file_LAMMPS_data ../../lammps_data_mixture_wrap.lmp \
  --file_LAMMPS_log log.lammps
mv VASP_2_ChIMES.xyzf LAMMPS.xyzf

