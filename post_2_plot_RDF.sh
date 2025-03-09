echo
echo "usage: \"do_script.sh  gr-1 gr-2 ... \""
echo "usage: \"do_script.sh  {all gr} \""
echo

if [ "$#" -eq 0 ]; then

    list_run=`ls -1vd gr_*.dat`

else

    for element in "$@" ; do
        list_run="${list_run} $element"
    done

fi


path="/p/lustre2/pham20/doing/TATB/99-MD-SCAN/MD-2500-3000K-1/1-8/"

for file in ${list_run} ; do
    echo $file
    xmgrace ~/tools/XMGRACE/tmp_RDF.agr  ${path}/$file $file
done
