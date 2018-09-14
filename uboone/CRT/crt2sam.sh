#!/bin/bash
for x in $(find $1 -type f -name CRTHits*.root);
do
  echo $x
  python ../extractor_dict_old.py ${x} > $(basename $1)_$(basename $(dirname ${x})).root.json
done
echo "====================="
echo 'created json files...'
echo "====================="

for x in $(find . -type f -name $(basename $1)_\*.root.json);
do
  python insert_art.py ${x} > ${x//crt/test_crt}
done
echo "======================================"
echo 'inserted family entry to json files...'
echo "======================================"
##=======================================================================
##for y in $(find . -type f -name test_crtdaq_20171130_18541\*.root.json);
##do
##var=$(grep -o 'CRTHits_ProdRun.*oot' $(basename ${y}))
##samweb retire-file $var
##done
###=======================================================================

for y in $(find . -type f -name test_$(basename $1)\*.root.json);
do
  echo $(basename ${y})
  samweb declare-file $(basename ${y})
  echo "declared CRTHits root file to SAM"
  offsprings=$(grep -o 'CRTHits_ProdRun.*oot' $(basename ${y}))
  echo $offsprings
  echo $1
  dir_address=$(dirname $(find $1 -type f -name $offsprings))
  echo "location of the offsprings $dir_address"
  samweb add-file-location $offsprings $dir_address
  echo "added location to SAM, should have gsiftp and xrootd URL now"
  echo "============================================================"
done
