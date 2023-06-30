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
cp -r $LARBATCH_DIR/python project_modules
touch project_modules/__init__.py

# Set experiment environment variables (not set by mrbsetenv, but needed by IFDH).

export EXPERIMENT=uboone
export SAM_EXPERIMENT=uboone

# This script runs the full mc+reco chain using standard released fcl files.

#input=$UBOONE_EXAMPLE_DATA_DIR/reco1/PhysicsRun-2018_10_11_9_59_56-0019527-00063_20181016T152442_ext_bnb_1_20190108T142430_optfilter_20190614T003502_reco1_postwcct_postdl_20190614T010503.root
input=$UBOONE_EXAMPLE_DATA_DIR/reco1/PhysicsRun-2018_10_11_9_59_56-0019527-00063_20181016T152442_ext_bnb_1_20190108T142430_optfilter_20190614T003502_reco1_postwcct_postdl_20190614T010503.root
fcl=test_sqlite_reco2_mcc9.fcl
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
