#! /bin/bash

# Skip for debug build.

if [[ x$MRB_QUALS =~ x.*debug.* ]]; then
  echo "Skipping for debug build."
  exit 0
fi
if [[ x`which lar` =~ x.*debug.* ]]; then
  echo "Skipping for debug build."
  exit 0
fi

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

# This script runs the full overlay chain using standard released fcl files.

input=$UBOONE_EXAMPLE_DATA_DIR/swizzled/PhysicsRun-2019_5_9_15_17_29-0022417-00341_20190528T225255_ext_bnb_8_20190528230111_merged.root
for fcl in standard_overlay_single_mu_driver.fcl wirecell_g4_uboone.fcl wirecell_detsim_overlay_uboone.fcl standard_overlay_uboone.fcl reco_uboone_mcc9_8_driver_overlay_stage1a.fcl reco_uboone_mcc9_8_driver_overlay_stage1b.fcl reco_uboone_mcc9_8_driver_overlay_stage1c.fcl wirecell_detsim_optical_overlay_uboone.fcl standard_overlay_optical_uboone.fcl reco_uboone_mcc9_8_driver_overlay_optical.fcl reco_uboone_mcc9_8_driver_overlay_stage2.fcl reco_uboone_data_mcc9_1_8_driver_poststage2_filters_Overlay.fcl
do
  output=`basename $fcl .fcl`.root
  out=`basename $fcl .fcl`.out
  err=`basename $fcl .fcl`.err
  cmd="lar --rethrow-all -c $fcl -s $input -o $output -n 1"
  echo $cmd
  $cmd > $out 2> $err
  stat=$?
  echo "Command finished with status $stat"
  if [ $stat -ne 0 ]; then
    exit $stat
  fi
  input=$output
done

# Done (success).

for root in *.root
do
  rootstat.py $root > ${root}.rootstat
  rm $root   # Clean up.
done

