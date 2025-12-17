E0=`grep NON VASP_51/VASP_2_ChIMES.xyzf | awk '{print $NF}'`

ids=`ls -d VASP_* | sed 's/VASP_//' | sort -n`
for id in $ids ; do
    fold="VASP_$id"
    E1=`grep NON VASP_$id/VASP_2_ChIMES.xyzf | awk '{print $NF}'`
    Ediff=`echo "scale=6; ($E1-1.0*$E0)/22.0" | bc -l`

    if (( $(echo "$Ediff > 1" | bc -l) )); then
        echo $fold
        echo $Ediff

        if [ -d "AAA" ]; then
            mv $fold AAA
        else
            mkdir AAA
            mv $fold AAA
        fi

    fi

done

