#! /bin/bash

# Make sure we can find config file in mrbsetenv environment.

if [ -d $PPFX_DIR/scripts ]; then
  export FW_SEARCH_PATH=$PPFX_DIR/scripts:$FW_SEARCH_PATH
fi

# Set up python path.

export PYTHONPATH=`pwd`:$UBUTIL_DIR/python:$LARBATCH_DIR/python:$PYTHONPATH
rm -rf project_modules
cp -r $LARBATCH_DIR/python project_modules
touch project_modules/__init__.py

# Save a copy of the environment (for debugging).

env > env.txt

# This script runs genie with flux type dk2nu and ppfx reweighting.

input=''
for fcl in test_genie_dk2nu.fcl ppfx_uboone.fcl
do
  output=`basename $fcl .fcl`.root
  out=`basename $fcl .fcl`.out
  err=`basename $fcl .fcl`.err
  if [ x$input = x ]; then
    cmd="lar --rethrow-all -c $fcl -o $output -n 2"
  else
    cmd="lar --rethrow-all -c $fcl -s $input -o $output -n 2"
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

# Done (success).

for root in *.root
do
  rootstat.py $root > ${root}.rootstat
  #rm $root   # Clean up.
done
