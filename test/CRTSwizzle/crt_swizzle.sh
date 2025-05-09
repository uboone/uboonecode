#! /bin/bash

# Save a copy of the environment (for debugging).

env > env.txt

# Exit if stash cache isn't mounted.

UBOONE_EXAMPLE_DATA_DIR=/pnfs/uboone/persistent/stash/uboone_example_data
if [ ! -d $UBOONE_EXAMPLE_DATA_DIR ]; then
  echo "Quittig because stash cache isn't available."
  exit
fi

# This script runs the standard swizzle fcl file.

input=$UBOONE_EXAMPLE_DATA_DIR/crtdaq/ProdRun20170603_161007-crt04.1.crtdaq.100k.part
fcl=CRTRaw2CRTHitsource1.fcl
out=`basename $fcl .fcl`.out
err=`basename $fcl .fcl`.err
cmd="lar --rethrow-all -c $fcl -s $input -n 1000000"
echo $cmd
$cmd > $out 2> $err
stat=$?
echo "Command finished with status $stat"
if [ $stat -ne 0 ]; then
  exit $stat
fi

# Done (success).

rm *.root   # Clean up.
