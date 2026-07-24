nimage=21

rm -rf neb.traj
python neb.py \
  --path_1 ../defect_energy/pu_vac_1/ \
  --path_2 ../defect_energy/pu_vac_2/ \
  --nimage $nimage
python neb_tools.py --nimage $nimage &> OUTPUT_1
mv energies.dat energies_1.dat

rm -rf neb.traj
python neb.py \
  --path_1 ../defect_energy/pu_oct/ \
  --path_2 ../defect_energy/pu_tet/ \
  --nimage $nimage
python neb_tools.py --nimage $nimage &> OUTPUT_2
mv energies.dat energies_2.dat
