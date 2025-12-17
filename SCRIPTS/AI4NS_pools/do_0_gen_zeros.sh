filePAIR="lammps_parameters_pairwise_llm105_buck.lmp"


molecule="ethyl_acetate"
fileoutput1="OUTPUT_${molecule}"
python gen_LJ_zero.py \
  --types_mol_1 1 2 3 4 5 6 \
  --types_mol_2 1 2 3 4 5 6 \
  &> $fileoutput1
./do_0_gen_pairwise.sh $fileoutput1 $filePAIR
mv NEW_PAIR ${filePAIR}_zero_$molecule
