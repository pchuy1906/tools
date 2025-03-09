module load python/2.7.18

file_xyz="input.xyzfes"
atom_types="C H N O"
polynomial_orders="12 8 3"
cutoff_distances="6 5 4"

python ~/tools/others/gen_fm_setup.py  \
--file_xyz ${file_xyz}  \
--atom_types ${atom_types}  \
--polynomial_orders ${polynomial_orders}  \
--cutoff_distances ${cutoff_distances} \
--delta_penalty 0.16 \
&> OUTPUT_0

