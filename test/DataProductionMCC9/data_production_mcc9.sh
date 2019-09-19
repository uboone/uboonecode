#! /bin/bash

# Save a copy of the environment (for debugging).

env > env.txt

# Exit if stash cache isn't mounted.

UBOONE_EXAMPLE_DATA_DIR=/cvmfs/uboone.osgstorage.org/stash/uboone_example_data
if [ ! -d $UBOONE_EXAMPLE_DATA_DIR ]; then
  echo "Quittig because stash cache isn't available."
  exit
fi

# Set wire cell path.

export WIRECELL_PATH=${UBOONEDATA_DIR}/WireCellData:${WIRECELL_FQ_DIR}/share/wirecell

# This script runs the full mc+reco chain using standard released fcl files.

input=$UBOONE_EXAMPLE_DATA_DIR/swizzled/PhysicsRun-2019_5_9_15_17_29-0022417-00341_20190528T225255_ext_bnb_8_20190528230111_merged.root
for fcl in reco_uboone_data_mcc9_8_driver_stage1.fcl reco_uboone_mcc9_8_driver_data_bnb_optical.fcl reco_uboone_data_mcc9_8_driver_stage2_reduced_beamOn.fcl standard_ana_uboone_data.fcl
do
  output=`basename $fcl .fcl`.root
  out=`basename $fcl .fcl`.out
  err=`basename $fcl .fcl`.err
  cmd="lar --rethrow-all -c $fcl -s $input -o $output -n 5"
  echo $cmd
  $cmd > $out 2> $err
  stat=$?
  echo "Command finished with status $stat"
  if [ $stat -ne 0 ]; then
    exit $stat
  fi
  input=$output
done

#for fcl in standard_larcv_uboone_data.fcl
#do
#  out=`basename $fcl .fcl`.out
#  err=`basename $fcl .fcl`.err
#  input=reco_uboone_data_mcc9_8_driver_stage2.root
#  cmd="lar --rethrow-all -c $fcl -s $input -n 5"
#  echo $cmd
#  $cmd > $out 2> $err
#  stat=$?
#  echo "Command finished with status $stat"
#  if [ $stat -ne 0 ]; then
#    exit $stat
#  fi
#done

# Done (success).

rm *.root   # Clean up.
