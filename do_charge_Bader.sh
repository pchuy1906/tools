# Bader
src="/p/lustre2/pham20/codes/bader_"
$src/chgsum.pl AECCAR0 AECCAR2 > OUT_Bader
$src/bader CHGCAR -ref CHGCAR_sum >> OUT_Bader
nions=`grep NIONS OUTCAR  | awk '{print $NF}'`
python ~/tools/VASP/export_ACF.py --natom $nions
