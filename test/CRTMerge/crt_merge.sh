#! /bin/bash

# Save a copy of the environment (for debugging).

env > env.txt

# This script runs the standard swizzle fcl file.

# Made up root file below needs to instead be a real file. -- EC, 17-Mar-2018
# We also require that CRTHits*.root files have been created for the time period of this TPCPMT file.
input=$UBOONE_EXAMPLE_DATA_DIR/crtdaq/PhysRunTPCPMT-swizzled-after_01-Dec-2017.root
fcl=testmerger.fcl
out=`basename $fcl .fcl`.out
err=`basename $fcl .fcl`.err
cmd="lar --rethrow-all -c $fcl -s $input -n 5"
echo $cmd
$cmd > $out 2> $err
stat=$?
echo "Command finished with status $stat"
if [ $stat -ne 0 ]; then
  exit $stat
fi
