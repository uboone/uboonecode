#! /bin/bash

# Save a copy of the environment (for debugging).

env > env.txt

# Set wire cell path.

export WIRECELL_PATH=${UBOONEDATA_DIR}/WireCellData:${WIRECELL_FQ_DIR}/share/wirecell

# This script runs the full mc+reco chain using standard released fcl files.

input=''
for fcl in prod_muminus_0.5-5.0GeV_25degf_uboone.fcl standard_g4_uboone.fcl standard_detsim_uboone.fcl reco_uboone_mcc8_driver_stage1_reduced.fcl reco_uboone_mcc8_driver_stage2_reduced.fcl standard_ana_uboone.fcl
do
  output=`basename $fcl .fcl`.root
  out=`basename $fcl .fcl`.out
  err=`basename $fcl .fcl`.err
  if [ x$input = x ]; then
    cmd="lar --rethrow-all -c $fcl -o $output -n 5"
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

for fcl in standard_larcv_uboone_mctruth.fcl standard_larcv_uboone.fcl
do
  out=`basename $fcl .fcl`.out
  err=`basename $fcl .fcl`.err
  if [ $fcl = standard_larcv_uboone_mctruth.fcl ]; then
    input=standard_detsim_uboone.root
  elif [ $fcl = standard_larcv_uboone.fcl ]; then
    input=reco_uboone_mcc8_driver_stage2_reduced.root
  fi
  cmd="lar --rethrow-all -c $fcl -s $input -n 5"
  echo $cmd
  $cmd > $out 2> $err
  stat=$?
  echo "Command finished with status $stat"
  if [ $stat -ne 0 ]; then
    exit $stat
  fi
done

