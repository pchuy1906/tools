echo
echo "usage: \"do_script.sh  file.xyz  cell_3/cell_9/NON_ORTHO \""
echo

if [ "$#" -ne 2 ]; then
  echo "unknown options. EXIT THE PROGRAM"
  exit 1
fi

filexyz=$1
cellopt=$2

echo "Input is: $file"

nat=`head -1 ${filexyz} | awk '{print $1}'`
echo
echo "The number of atom is " $nat

cat > qe.in << EOF
&control
  calculation='scf'
  prefix='water_SiO2'
  pseudo_dir='/p/lustre2/pham20/SiO2/QE_from_Marcos/'
  outdir='./pwscf'
  restart_mode='from_scratch'
  nstep=20000
  tprnfor=.true.
  tstress=.true.
/
&system
  ibrav=0
  nat=$nat
  ntyp=3,
  ecutwfc=110
  input_dft='SCAN'
occupations='smearing'
smearing='fd'
degauss=0.002
/
&electrons
  electron_maxstep = 1000
  mixing_beta = 0.2
  conv_thr = 1.D-6
/
&ions
/
&cell
/
ATOMIC_SPECIES
O    1.0  O.tm.pbe.UPF 
H    1.0  H.tm.pbe.UPF 
Si   1.0  Si.pbe-n-nc.UPF
EOF

echo "CELL_PARAMETERS {angstrom}" >> qe.in

if [ "$cellopt" == "cell_9" ]; then
    head -2 ${filexyz} | tail -1 | awk '{for(i=1;i<=NF;i+=3) printf "%15.9f %15.9f %15.9f\n", $i, $(i+1), $(i+2)}' >> qe.in
elif [ "$cellopt" == "cell_3" ]; then
    a=`head -2 ${filexyz} | tail -1 | awk '{print $1}'`
    b=`head -2 ${filexyz} | tail -1 | awk '{print $2}'`
    c=`head -2 ${filexyz} | tail -1 | awk '{print $3}'`
    echo $a 0.0 0.0 >> qe.in
    echo 0.0 $b 0.0 >> qe.in
    echo 0.0 0.0 $c >> qe.in
elif [ "$cellopt" == "NON_ORTHO" ]; then
    head -2 ${filexyz} | tail -1 | sed 's/NON_ORTHO//g' | awk '{for(i=1;i<=NF;i+=3) printf "%15.9f %15.9f %15.9f\n", $i, $(i+1), $(i+2)}' >> qe.in
fi

echo "" >> qe.in
echo "ATOMIC_POSITIONS {angstrom}" >> qe.in
tail -$nat ${filexyz} | awk '{printf "%5s %15.9f %15.9f %15.9f\n", $1, $2, $3, $4}' >> qe.in
