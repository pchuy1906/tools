cp2k="/g/g92/pham20/codes/cp2k/cp2k-7.1/"

for element in C O N F S Zn Yb ; do
    echo $element
    grep -e "$element " ${cp2k}/data/BASIS_MOLOPT* | grep -e "DZVP"
done
