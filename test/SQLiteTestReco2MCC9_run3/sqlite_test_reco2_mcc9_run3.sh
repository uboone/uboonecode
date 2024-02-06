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

export PYTHONPATH=`pwd`:$UBUTIL_DIR/python:$LARBATCH_DIR/python:$PYTHONPATH
rm -rf project_modules
if [ -d $LARBATCH_DIR/python/project_modules ]; then
  cp -r $LARBATCH_DIR/python/project_modules project_modules
else
  cp -r $LARBATCH_DIR/python project_modules
fi
touch project_modules/__init__.py

# Set experiment environment variables (not set by mrbsetenv, but needed by IFDH).

export EXPERIMENT=uboone
export SAM_EXPERIMENT=uboone

# This script runs the full mc+reco chain using standard released fcl files.

input=$UBOONE_EXAMPLE_DATA_DIR/reco1/PhysicsRun-2018_7_2_21_23_1-0017505-00141_20180714T132533_ext_bnb_3_20181207T001756_optfilter_20181224T035203_reco1_postwcct_postdl_20181224T035451.root
fcl=test_sqlite_reco2_mcc9_run3.fcl
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

# Done (success).

for root in *.root
do
  rootstat.py $root > ${root}.rootstat
  #rm $root   # Clean up.
done
