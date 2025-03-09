for ID in  `cat LIST_STRUCTURE` ; do
  echo $ID

  cd VASP_${ID}
    /g/g92/pham20/tools/take_xyzfes.sh
  cd ..

  echo
  echo
done
