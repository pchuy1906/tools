for alpha in 10000 5000 1000 500 100 50 10 5 1 0.5 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 ; do

    #rm -rf alpha_${alpha}; mkdir alpha_${alpha}
    #python /g/g92/pham20/codes/quartz/fmatch-2018-r456/generalized-branch-4b/src/lsq2.py  --algorithm lassolars --alpha ${alpha}  --weights  weights.dat    > alpha_${alpha}/params.txt_1e-2_LASSO_LARS 
    #mv force.txt  Xdata  alpha_${alpha}

    #rm -rf alpha_${alpha}; mkdir alpha_${alpha}
    #cd alpha_${alpha}
    #    cp ../head_job2.csh  job2.sh
    #    echo "python /g/g92/pham20/codes/quartz/fmatch-2018-r456/generalized-branch-4b/src/lsq2.py  --A  ../A.txt  --b ../b.txt --algorithm lassolars --alpha ${alpha}  --weights  ../weights.dat    > params.txt_1e-2_LASSO_LARS"  >> job2.sh
    #    msub  job2.sh
    #cd ..

    max_abs=`python ~/tools/max_abs.py   --file_input  alpha_${alpha}/Xdata | awk '{print $1}'`
    norm_LA=`python ~/tools/max_abs.py   --file_input  alpha_${alpha}/Xdata | awk '{print $2}'`
    num_var=`grep -e "number of fitting vars"  alpha_${alpha}/params.txt_1e-2_LASSO_LARS  | awk '{print $NF}'`
    paste DFT.FORCE_ENERGY_STRESS  alpha_${alpha}/force.txt | grep -e "FORCE" | awk '{print $1, $3}' > COMPARE_F.dat
    RMS_F=`python ~/tools/RMS.py  --file_input  COMPARE_F.dat | tail -1 | awk '{print $NF}'`
    echo ${max_abs} ${RMS_F} ${num_var} ${norm_LA} ${alpha}

done
