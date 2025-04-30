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

UBOONE_EXAMPLE_DATA_DIR=/pnfs/uboone/persistent/stash/uboone_example_data
if [ ! -d $UBOONE_EXAMPLE_DATA_DIR ]; then
  echo "Quittig because stash cache isn't available."
  exit
fi

# Set wire cell path.

export WIRECELL_PATH=${UBOONEDATA_DIR}/WireCellData:${WIRECELL_FQ_DIR}/share/wirecell

# Set up python path.

export PYTHONPATH=`pwd`:$UBUTIL_DIR/python:$LARBATCH_DIR/python:$PYTHONPATH
rm -rf project_modules
cp -r $LARBATCH_DIR/python project_modules
touch project_modules/__init__.py

# Set experiment environment variables (not set by mrbsetenv, but needed by IFDH).

export EXPERIMENT=uboone
export SAM_EXPERIMENT=uboone

# This script runs the nugraph chain using standard released fcl files.

input=$UBOONE_EXAMPLE_DATA_DIR/mcreco1/PhysicsRun-2020_2_24_2_35_6-0025487-00006_20200419T164054_ext_unbiased_1_20200422T181434_merged_gen_fhc_20241110T063112_g4_detsim_mix_r1a_r1b_20241110T065513_reco1c_20241111T074036_merged.root
for fcl in run_reboone.fcl wirecell_reg4_LightPropTime_LY.fcl wirecell_detsim_optical_overlay_uboone.fcl standard_overlay_notpc_uboone.fcl reco_uboone_mcc9_8_driver_overlay_optical.fcl reco_uboone_mcc9_10_driver_overlay_stage2.fcl run_fullshowerreco.fcl
do
  output=`basename $fcl .fcl`.root
  out=`basename $fcl .fcl`.out
  err=`basename $fcl .fcl`.err
  cmd="lar --rethrow-all -c $fcl -s $input -o $output -n 2"
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
