echo "usage: \"do_script.sh    fold1_ambient    fold2_higher    fold3_lower \""
if [ "$#" -ne 3 ]; then
    echo "unknown options. EXIT THE PROGRAM"
    exit 1
fi

fold1=$1
fold2=$2
fold3=$3

curr_dir=`pwd`
for fold in $fold1 $fold2 $fold3 ; do
    cd $fold
        ~/tools/LAMMPS/scripts/do_collect_LAMMPS.sh  &> OUTPUT_1

        head -n -1 AAA_Volume.dat   > tmp_AAA_Volume.dat
        head -n -1 AAA_Temp.dat     > tmp_AAA_Temp.dat
        head -n -1 AAA_Press.dat    > tmp_AAA_Press.dat
        head -n -1 AAA_TotEng.dat   > tmp_AAA_TotEng.dat
        ndata=`wc tmp_AAA_Press.dat | awk '{print $1}'`
        nhalf=$(($ndata/2))
        tail -${nhalf} tmp_AAA_Volume.dat   > AAA2_Volume.dat
        tail -${nhalf} tmp_AAA_Temp.dat     > AAA2_Temp.dat
        tail -${nhalf} tmp_AAA_Press.dat    > AAA2_Press.dat
        tail -${nhalf} tmp_AAA_TotEng.dat   > AAA2_TotEng.dat


        python /g/g92/pham20/tools/others/block_ave.py --file_input AAA2_Volume.dat --nblock 8 > OUTPUT_2
        grep -e "average" _OOO | awk '{print $(NF-1), $NF}' >  average_data.dat
    
        python /g/g92/pham20/tools/others/block_ave.py --file_input AAA2_Temp.dat   --nblock 8 > OUTPUT_2
        grep -e "average" _OOO | awk '{print $(NF-1), $NF}' >> average_data.dat
    
        python /g/g92/pham20/tools/others/block_ave.py --file_input AAA2_Press.dat  --nblock 8 > OUTPUT_2
        grep -e "average" _OOO | awk '{print $(NF-1), $NF}' >> average_data.dat
    
        python /g/g92/pham20/tools/others/block_ave.py --file_input AAA2_TotEng.dat --nblock 8 > OUTPUT_2
        grep -e "average" _OOO | awk '{print $(NF-1), $NF}' >> average_data.dat

    cd ${curr_dir}
done

python ~/tools/others/calculate_Hugo_full_2.py \
--file0 $fold1/average_data.dat  \
--file1 $fold2/average_data.dat  \
--file2 $fold3/average_data.dat  \
--unit_P bars \
--unit_E eV

