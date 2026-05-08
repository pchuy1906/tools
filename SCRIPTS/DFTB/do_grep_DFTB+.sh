echo "usage: \"do_script.sh    DFTB_output  0 0\""
if [ "$#" -ne 3 ]; then
    echo "unknown options. EXIT THE PROGRAM"
    exit 1
fi

fileout=$1
n0=$2
nskip=$3


#Total MD Energy:                  -327.4859046789 H        -8911.3449 eV
grep -e "Total MD Energy:" ${fileout}  | awk '{print $(NF-1)}' >  TotEng.dat 
#Total Energy:                     -327.7635734670 H        -8918.9006 eV
#grep -e "Total Energy:"    ${fileout}  | awk '{print $(NF-1)}' >  TotEng.dat
#Force related energy:             -493.4950894662 H       -13428.6846 eV
#grep -e "Force related energy:"    ${fileout}  | awk '{print $(NF-1)}' >  TotEng.dat



#MD Temperature:                      0.0008846317 H          279.3443 K
grep -e "MD Temperature"   ${fileout}  | awk '{print $(NF-1)}' >  Temp.dat

#Pressure:                           -0.599780E-04 au    -0.176461E+10 Pa
grep -e "Pressure"         ${fileout}  | awk '{print $(NF-1)}' >  Press.dat
awk '{print $1*0.000000001}' Press.dat > TMP ; mv TMP Press.dat

#Volume:                              0.813383E+05 au^3   0.120531E+05 A^3
grep -e "Volume"           ${fileout}  | awk '{print $(NF-1)}' >  Volume.dat

for file in TotEng.dat Temp.dat Press.dat Volume.dat ; do

    awk '{print NR+'"$n0"' " " $s}'   $file  > TMP

    cp TMP AAA2_${file}

    if [ ${nskip} -ge 3 ]; then
        tail -n +"${nskip}" TMP > AAA2_${file}
    fi

    #rm -rf $file

done

rm -rf TMP

