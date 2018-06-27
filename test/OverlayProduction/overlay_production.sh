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

# This script runs the full overlay chain using standard released fcl files.

input=$UBOONE_EXAMPLE_DATA_DIR/swizzled/PhysicsRun-2016_3_14_9_22_21-0005432-00021_20160322T065603_ext_bnb_20160323T041757_merged.root
for fcl in standard_overlay_single_mu_driver.fcl standard_g4_overlay_uboone.fcl standard_detsim_nonoise_overlay_uboone.fcl standard_overlay_uboone.fcl reco_uboone_mcc8_driver_overlay_stage1a.fcl reco_uboone_mcc8_driver_overlay_stage1b.fcl reco_uboone_mcc8_driver_overlay_stage1c.fcl reco_uboone_mcc8_driver_overlay_stage2.fcl
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
