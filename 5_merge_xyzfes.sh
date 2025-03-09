rm -rf merged_input.xyzfes  merged_DFT.FORCE_ENERGY_STRESS
for ID in  `cat LIST_STRUCTURE` ; do
    cat VASP_${ID}/input.xyzfes >> merged_input.xyzfes
    cat VASP_${ID}/DFT.FORCE_ENERGY_STRESS  >> merged_DFT.FORCE_ENERGY_STRESS
done
