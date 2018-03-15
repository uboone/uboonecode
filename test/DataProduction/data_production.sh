#! /bin/bash

# Save a copy of the environment (for debugging).

env > env.txt

# This script runs the full mc+reco chain using standard released fcl files.

input=$UBOONE_EXAMPLE_DATA_DIR/swizzled/PhysicsRun-2016_3_14_9_22_21-0005432-00021_20160322T065603_ext_bnb_20160323T041757_merged.root
for fcl in reco_uboone_data_mcc8_driver_stage1_reduced.fcl reco_uboone_data_mcc8_driver_stage2_reduced.fcl standard_ana_uboone_data.fcl
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

for fcl in standard_larcv_uboone_data.fcl
do
  out=`basename $fcl .fcl`.out
  err=`basename $fcl .fcl`.err
  input=reco_uboone_data_mcc8_driver_stage2_reduced.root
  cmd="lar --rethrow-all -c $fcl -s $input -n 5"
  echo $cmd
  $cmd > $out 2> $err
  stat=$?
  echo "Command finished with status $stat"
  if [ $stat -ne 0 ]; then
    exit $stat
  fi
done
