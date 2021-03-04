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

# Set up python path.

export PYTHONPATH=$UBUTIL_DIR/python:$PYTHONPATH

# Set experiment environment variables (not set by mrbsetenv, but needed by IFDH).

export EXPERIMENT=uboone
export SAM_EXPERIMENT=uboone

# This script runs the full mc+reco chain using standard released fcl files.

input=$UBOONE_EXAMPLE_DATA_DIR/swizzled/PhysicsRun-2019_5_9_15_17_29-0022417-00198_20190620T162349_bnb_2_20190620122031_merged.root
for fcl in run_merge_beamdata.fcl reco_uboone_data_mcc9_8_driver_stage1.fcl reco_uboone_mcc9_8_driver_data_bnb_optical.fcl reco_uboone_data_mcc9_8_driver_stage2_reduced_beamOn.fcl reco_uboone_data_mcc9_1_8_driver_poststage2_filters_beamOn.fcl
do
  output=`basename $fcl .fcl`.root
  out=`basename $fcl .fcl`.out
  err=`basename $fcl .fcl`.err
  if [ $fcl = reco_uboone_data_mcc9_1_8_driver_poststage2_filters_beamOn.fcl ]; then
    cmd="lar --rethrow-all -c $fcl -s $input -n 5"
  else
    cmd="lar --rethrow-all -c $fcl -s $input -o $output -n 5"
  fi
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

for root in *.root
do
  rootstat.py $root > ${root}.rootstat
  rm $root   # Clean up.
done
